# Summary of one state's pair-check units (pc/$PC_LABEL/units): counts, every multi-root cluster
# in pc/$PC_LABEL/clusters.csv, every root absent from production in pc/$PC_LABEL/absent.csv.
# Tolerances: the mesh study's observation tolerances (0.15 dB, 1.5°, 0.019 in s2 and s3).
using Pkg
const ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(ROOT)
using JLD2, Printf

const TOL = (amp=0.15, phase=1.5, s2=0.019, s3=0.019)

function main()
    label = ENV["PC_LABEL"]
    expected = parse(Int, get(ENV, "PC_EXPECTED", "0"))
    dir = joinpath(ROOT, "pc", label)
    files = filter(f -> endswith(f, ".jld2"), readdir(joinpath(dir, "units"); join=true))
    classes = Dict{String,Int}()
    nseg = nraw = nabsent = nunmatched = nexceed = 0
    worst = (amp=0.0, phase=0.0, s2=0.0, s3=0.0)
    open(joinpath(dir, "clusters.csv"), "w") do cio
    open(joinpath(dir, "absent.csv"), "w") do aio
        println(cio, "member,path,segment,class,n_raw,boxes,min_gap_rad,count_radius_rad,zeros,located,verified,verify_radius_rad,winding,reliable,poles_inside,poles_near,re_deg,im_deg")
        println(aio, "member,path,segment,re_deg,im_deg,cluster_class,d_amp_dB,d_phase_deg,d_s2,d_s3")
        for f in files
            u = load(f)
            r, m, path = u["result"], u["member"], u["path"]
            Δ = r.fields === missing ? missing : r.fields.delta
            if Δ !== missing
                worst = map((w, d) -> max(w, abs(d)), worst, Δ)
                any(abs(getfield(Δ, k)) > getfield(TOL, k) for k in keys(TOL)) && (nexceed += 1)
            end
            for (k, s) in enumerate(r.segments)
                nseg += 1; nraw += length(s.raw); nunmatched += length(s.unmatched)
                for c in s.clusters
                    classes[c.class] = get(classes, c.class, 0) + 1
                    @printf(cio, "%d,%s,%d,%s,%d,%s,%.3e,%.3e,%d,%d,%s,%s,%.4f,%s,%d,%d,%.6f,%.6f\n", m, path, k, c.class,
                        length(c.raw), join(c.boxes, ";"), minimum(c.gaps), c.radius, c.zeros, length(c.roots),
                        join(c.verified, ";"), join((@sprintf("%.3e", r) for r in c.vradius), ";"),
                        c.winding, c.reliable, c.poles_inside, c.poles_near_circle,
                        rad2deg(real(c.center)), rad2deg(imag(c.center)))
                end
                for θ in s.absent
                    nabsent += 1
                    cl = something(findfirst(c -> any(d -> abs(d - θ) < 1e-7, c.roots), s.clusters), 0)
                    @printf(aio, "%d,%s,%d,%.6f,%.6f,%s,%s\n", m, path, k, rad2deg(real(θ)), rad2deg(imag(θ)),
                        cl == 0 ? "none" : s.clusters[cl].class,
                        Δ === missing ? ",,," : join((@sprintf("%.4g", getfield(Δ, q)) for q in keys(TOL)), ","))
                end
            end
        end
    end
    end
    @printf("%s: %d units%s, %d segments, %d raw roots\n", label, length(files),
        expected > 0 ? " of $expected" : "", nseg, nraw)
    println("clusters (raw roots within the radius): ", isempty(classes) ? "none" :
        join(("$k $v" for (k, v) in sort!(collect(classes))), ", "))
    @printf("roots absent from production: %d; production roots not found by the check: %d\n", nabsent, nunmatched)
    @printf("units with absent roots outside tolerance: %d; largest |Δ|: %.3g dB, %.3g°, s2 %.3g, s3 %.3g\n",
        nexceed, worst.amp, worst.phase, worst.s2, worst.s3)
end

main()
