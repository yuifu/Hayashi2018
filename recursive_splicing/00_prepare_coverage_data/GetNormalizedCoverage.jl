using Bio.Align
using Glob
using GZip

function coveragemap(reader, chrom, region)
    cov = zeros(Int, length(region))

    aln = Bio.Align.Alignment("")
    
    for rec in intersect(reader, chrom, region)
        if !ismapped(rec)
            continue
        end

        aln = alignment(rec)
        for i in 1:seqlength(rec)
            j, op = seq2ref(aln, i)
            if ismatchop(op) && j in region
                cov[j - first(region) + 1] += 1
            end
        end
    end
    return cov
end


function batchOpenBamReader(bamFiles)
    map(bf -> open(BAMReader, bf, index=bf * ".bai"), bamFiles)
end

function batchCloseBamReader(bamReaders)
    map(x -> close(x), bamReaders)
end

# function getCoverage(bamReaders, chrom, region)
#     mat = zeros(length(bamReaders), length(region))
#     for i in 1:length(bamReaders)
#         mat[i,:] = coveragemap(bamReaders[i], chrom, region)
#     end
#     return mat
# end

function getNormalizedCoverageSum(bamReaders, chrom, region, normFac)
    arr = zeros(length(region))
    cov = zeros(Int, length(region))
    for i in 1:length(bamReaders)
        cov = coveragemap(bamReaders[i], chrom, region)
        for j in 1:length(cov)
            arr[j] +=  cov[j] / normFac[i]
        end
    end
    return arr
end

function getNormalizedCoverageBase(bamReaders, chrom, region, normFac)
    mat = zeros(length(bamReaders), length(region))
    for i in 1:length(bamReaders)
        mat[i,:] = coveragemap(bamReaders[i], chrom, region) / normFac[i]
    end
    return vec(sum(mat, 1))
end


function parsePathBmFiles(pathBamFiles)
    bamFiles = Array{Any,1}()
    normFac = Array{Float64,1}()

    open(pathBamFiles) do f
        readline(f)
        for line in eachline(f)
            arr = split(chomp(line), "\t")
            push!(bamFiles, arr[1])
            push!(normFac, parse(Float64, arr[2]))
        end
    end

    return bamFiles, normFac
end


function parsePathRegions(pathRegions)
    I_CHROM = 1
    I_LPOS = 2
    I_RPOS = 3
    I_whichIsExonAnchor = 5
    I_NAME = 4
    I_REF_LPOS = 8
    I_REF_RPOS = 9

    ids = Array{Any,1}()
    chroms = Array{Any,1}()
    regions = Array{UnitRange{Int64},1}()
    rsss = Array{Int64,1}()

    rss = 0

    open(pathRegions) do f
        readline(f)
        for line in eachline(f)
            arr = split(chomp(line), "\t")

            if arr[I_whichIsExonAnchor] == "1"
                rss = parse(Int,arr[I_RPOS])
            else
                rss = parse(Int,arr[I_LPOS]) + 1
            end

            region = (parse(Int,arr[I_REF_LPOS])+1):(parse(Int,arr[I_REF_RPOS]))
            id = arr[I_NAME]

            push!(ids, id)
            push!(regions, region)
            push!(rsss, rss)
            push!(chroms, arr[I_CHROM])
        end
    end

    return ids, chroms, regions, rsss
end


function writeResult(ofile, region, cov, rss)
    RSS = ""
    GZip.open(ofile, "w") do fw
        write(fw, "bp\tdepth\tRSS\n")
        for i in 1:length(region)
            if region[i] < rss
                RSS = "TRUE"
            else
                RSS = "FALSE"
            end

            x = @sprintf "%d\t%.5f\t%s\n" region[i] cov[i] RSS
            write(fw, x)
        end
    end
end


function getNormalizedCoverage(pathBamFiles, pathRegions, odirPrefix)

    bamFiles, normFac = parsePathBmFiles(pathBamFiles)
    println(bamFiles)
    ids, chroms, regions, rsss = parsePathRegions(pathRegions)

    bamReaders = batchOpenBamReader(bamFiles)

    for i in 1:length(ids)
        println((@sprintf "Processing %s" ids[i]))

        id = ids[i]
        chrom = chroms[i]
        region = regions[i]
        rss = rsss[i]

        cov = getNormalizedCoverageSum(bamReaders, chrom, region, normFac)

        ofile = @sprintf "%s/%s.cov.txt.gz" odirPrefix id

        writeResult(ofile, region, cov, rss)
    end

    batchCloseBamReader(bamReaders)
end





