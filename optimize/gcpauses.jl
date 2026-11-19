using Printf, Statistics

const PAUSE = r"GC: pause ([0-9.]+)ms\. collected ([0-9.]+)MB\. (full|incr) ?(recollect)?"
const MARK = r"== GCLOG (start|end) ([^\r\n]+)"

pct(v, q) = isempty(v) ? NaN : quantile(v, q)

function describe(label, ms, mb, full, recollect)
    f, i = ms[full], ms[.!full]
    @printf("%-12s pauses %6d (full %5d, recollect %4d); total %8.1f s (full %7.1f s); collected %8.1f GiB\n",
        label, length(ms), count(full), count(recollect), sum(ms; init=0.0)/1e3, sum(f; init=0.0)/1e3,
        sum(mb; init=0.0)/1024)
    for (name, v) in (("incr", i), ("full", f))
        @printf("    %-4s ms: median %8.1f  p90 %8.1f  p99 %8.1f  max %8.1f\n",
            name, pct(v, 0.5), pct(v, 0.9), pct(v, 0.99), isempty(v) ? NaN : maximum(v))
    end
end

"""
    gcpauses(logfile)

Summarize `GC.enable_logging(true)` output in `logfile`, per `== GCLOG start/end <label>` window,
or for the whole file if there are no markers.
"""
function gcpauses(logfile)
    text = read(logfile, String)
    marks = [(m.offset, m.captures[1], strip(m.captures[2])) for m in eachmatch(MARK, text)]
    pauses = [(m.offset, parse(Float64, m.captures[1]), parse(Float64, m.captures[2]),
        m.captures[3] == "full", m.captures[4] !== nothing) for m in eachmatch(PAUSE, text)]
    windows = Tuple{String,Int,Int}[]
    for (o, kind, label) in marks
        kind == "start" && push!(windows, (label, o, typemax(Int)))
        if kind == "end"
            k = findlast(w -> w[1] == label && w[3] == typemax(Int), windows)
            k === nothing || (windows[k] = (label, windows[k][2], o))
        end
    end
    isempty(windows) && push!(windows, ("whole log", 0, typemax(Int)))
    for (label, a, b) in windows
        sel = [p for p in pauses if a < p[1] < b]
        describe(label, [p[2] for p in sel], [p[3] for p in sel], Bool[p[4] for p in sel], Bool[p[5] for p in sel])
    end
end

gcpauses(only(ARGS))
