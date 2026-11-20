# Runs Main.jl as configured by the environment and saves timing, GC and provenance to
# $NEW_FOLDER/g5run.jld2. Workers started by Main.jl (DISTRIBUTED=true) report their GC since start.
using Pkg
const ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(ROOT)
using Distributed, JLD2, Logging, Dates, Printf

"""Forwards every message to `inner` and records the `Δt` of Main.jl's per-iteration "Elapsed" message."""
struct ElapsedRecorder <: AbstractLogger
    inner::AbstractLogger
    seconds::Vector{Float64}
end
Logging.min_enabled_level(l::ElapsedRecorder) = Logging.min_enabled_level(l.inner)
Logging.shouldlog(l::ElapsedRecorder, args...) = Logging.shouldlog(l.inner, args...)
Logging.catch_exceptions(l::ElapsedRecorder) = Logging.catch_exceptions(l.inner)
function Logging.handle_message(l::ElapsedRecorder, level, message, args...; kwargs...)
    message == "Elapsed" && haskey(kwargs, :Δt) && push!(l.seconds, Dates.toms(kwargs[:Δt]) / 1000)
    Logging.handle_message(l.inner, level, message, args...; kwargs...)
end

const RECORDER = ElapsedRecorder(global_logger(), Float64[])
global_logger(RECORDER)

const CONFIG = (; (Symbol(k) => get(ENV, k, "") for k in ("PATHS", "TOD", "ENS_SIZE", "ITRS",
    "NEW_FOLDER", "DISTRIBUTED", "NPROCS", "PROC_THREADS", "PROC_FLAGS"))...)
const PROV = (; commit=readchomp(`git -C $ROOT rev-parse HEAD`),
    dirty=!isempty(readchomp(`git -C $ROOT status --porcelain --untracked-files=no`)),
    host=gethostname(), nthreads=Threads.nthreads(), julia=string(VERSION), started=now())
@info "g5run" CONFIG PROV

const RUN = @timed include(joinpath(ROOT, "Main.jl"))

const WORKERS = nprocs() == 1 ? NamedTuple[] : [begin
        n = remotecall_fetch(Base.gc_num, w)
        (; pid=w, nthreads=remotecall_fetch(Threads.nthreads, w), gctime=n.total_time / 1e9,
            bytes=Base.gc_total_bytes(n), pauses=Int(n.pause), full=Int(n.full_sweep),
            maxrss=remotecall_fetch(Sys.maxrss, w))
    end for w in workers()]

const RESULT = (; wall=RUN.time, iterations=copy(RECORDER.seconds), gctime=RUN.gctime,
    bytes=RUN.bytes, pauses=Int(RUN.gcstats.pause), full=Int(RUN.gcstats.full_sweep),
    maxrss=Sys.maxrss(), workers=WORKERS)

@printf("== total %.1f s; iterations %s s\n", RESULT.wall, join(round.(RESULT.iterations; digits=1), ", "))
@printf("   main process: GC %.1f s (%.1f%%), %.1f GiB allocated, pauses %d, full %d, maxrss %.1f GiB\n",
    RESULT.gctime, 100RESULT.gctime / RESULT.wall, RESULT.bytes / 2^30, RESULT.pauses, RESULT.full,
    RESULT.maxrss / 2^30)
for w in WORKERS
    @printf("   worker %d (%d threads): GC %.1f s (%.1f%% of total), %.1f GiB allocated, pauses %d, full %d, maxrss %.1f GiB\n",
        w.pid, w.nthreads, w.gctime, 100w.gctime / RESULT.wall, w.bytes / 2^30, w.pauses, w.full, w.maxrss / 2^30)
end
jldsave(joinpath(ROOT, CONFIG.NEW_FOLDER, "g5run.jld2"); config=CONFIG, provenance=PROV, result=RESULT)
