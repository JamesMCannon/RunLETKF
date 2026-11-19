using Logging

"""
    PolishMesh(inner; maxit=8, xtol=1e-8, maxstep=5e-4, retry_tolerance=1e-9, retry_maxit=20,
               filter_radius=3e-4, filter_points=64, guard90=1e-3, dedupe=1e-7)

Roots of `inner`, finalized as `AVIDModes.finalize`: each is Newton-polished (with a
tight-integration retry); one that does not converge is kept unpolished only if the winding number
of the modal function on a `filter_radius` circle around it is +1; polished roots within `guard90`
rad of π/2 and duplicates within `dedupe` rad are dropped. Every root is logged to `EVENTS` and
every search to `SEARCHES`, tagged with the task's `:p1` context `(member, seq)`.
"""
Base.@kwdef struct PolishMesh{M}
    inner::M
    maxit::Int = 8
    xtol::Float64 = 1e-8
    maxstep::Float64 = 5e-4
    retry_tolerance::Float64 = 1e-9
    retry_maxit::Int = 20
    filter_radius::Float64 = 3e-4
    filter_points::Int = 64
    guard90::Float64 = 1e-3
    dedupe::Float64 = 1e-7
end
PolishMesh(inner; kw...) = PolishMesh(; inner, kw...)

const EventT = NamedTuple{(:member, :seq, :raw, :polished, :outcome, :winding, :abort, :loss),
    Tuple{Int,Int,ComplexF64,ComplexF64,Symbol,Float64,Symbol,Float64}}
const SearchT = NamedTuple{(:member, :seq, :nraw, :nkept, :abort, :t_search, :t_polish, :bytes_polish),
    Tuple{Int,Int,Int,Int,Symbol,Float64,Float64,Int}}
const EVENTS = EventT[]
const SEARCHES = SearchT[]
const EVENT_LOCK = ReentrantLock()

"""Logger that records a RootsAndPoles `maxiters`/`maxnodes` abort and forwards every message."""
struct AbortLogger <: AbstractLogger
    parent::AbstractLogger
    abort::Base.RefValue{Symbol}
end
Logging.min_enabled_level(l::AbortLogger) = Logging.min_enabled_level(l.parent)
Logging.shouldlog(::AbortLogger, args...) = true
Logging.catch_exceptions(l::AbortLogger) = Logging.catch_exceptions(l.parent)
function Logging.handle_message(l::AbortLogger, level, message, args...; kwargs...)
    s = string(message)
    occursin("`maxiters` reached", s) && (l.abort[] = :maxiters)
    occursin("`maxnodes` reached", s) && (l.abort[] = :maxnodes)
    return Logging.handle_message(l.parent, level, message, args...; kwargs...)
end

loss_dBMm(θ, k) = -20/log(10) * k * imag(sin(θ)) * 1e6

gcbytes() = (b = Ref{Int64}(0); Base.gc_bytes(b); b[])

function polish(me, θ0, mesh::PolishMesh; params, maxit)
    θ = θ0
    try
        for _ in 1:maxit
            dF, R, Rg = LMP.solvedmodalequation(θ, me; params)
            δ = LMP.modalequation(R, Rg) / dF
            isfinite(δ) || return θ0, false
            θ -= δ
            abs(θ - θ0) > mesh.maxstep && return θ0, false
            abs(δ) < mesh.xtol && return θ, true
        end
    catch
    end
    return θ0, false
end

function winding(me, θ0, r; n, params)
    z = θ0 .+ r .* cis.(range(0, 2π; length=n+1))
    f = [LMP.solvemodalequation(θ, me; params) for θ in z]
    return sum(angle.(f[2:end] ./ f[1:end-1])) / 2π
end

function LMP.findmodes(me::LMP.ModeEquation, mesh::PolishMesh; params=LMPParams())
    ctx = get(task_local_storage(), :p1, (member=0, seq=Ref(0)))
    ctx.seq[] += 1
    abort = Ref(:none)
    t0 = time_ns()
    raw = with_logger(AbortLogger(current_logger(), abort)) do
        ComplexF64[ea.θ for ea in LMP.findmodes(me, mesh.inner; params)]
    end
    t1, b1 = time_ns(), gcbytes()
    tight = LMPParams(params; integrationparams=LMP.IntegrationParams(solver=LMP.Vern7(),
        tolerance=mesh.retry_tolerance))
    kept, evs = ComplexF64[], EventT[]
    for θ0 in raw
        θ, ok = polish(me, θ0, mesh; params, maxit=mesh.maxit)
        outcome = :converged
        if !ok
            θ, ok = polish(me, θ0, mesh; params=tight, maxit=mesh.retry_maxit)
            outcome = :retry_converged
        end
        w = NaN
        if !ok
            w = winding(me, θ0, mesh.filter_radius; n=mesh.filter_points, params)
            outcome = isfinite(w) && round(Int, w) == 1 ? :unconverged_kept : :filter_drop
        end
        if outcome != :filter_drop
            if abs(θ - π/2) < mesh.guard90
                outcome = :guard_drop
            elseif any(x -> abs(x - θ) < mesh.dedupe, kept)
                outcome = :dup_drop
            else
                push!(kept, θ)
            end
        end
        push!(evs, (; member=ctx.member, seq=ctx.seq[], raw=θ0, polished=θ, outcome, winding=w,
            abort=abort[], loss=loss_dBMm(θ0, me.frequency.k)))
    end
    t2, b2 = time_ns(), gcbytes()
    lock(EVENT_LOCK) do
        append!(EVENTS, evs)
        push!(SEARCHES, (; member=ctx.member, seq=ctx.seq[], nraw=length(raw), nkept=length(kept),
            abort=abort[], t_search=(t1 - t0)/1e9, t_polish=(t2 - t1)/1e9, bytes_polish=Int(b2 - b1)))
    end
    return LMP.EigenAngle.(kept)
end
