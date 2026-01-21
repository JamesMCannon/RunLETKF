
function circular_mean(x; period=4)
    θ = 2π .* x ./ period
    μθ = atan(mean(sin.(θ)), mean(cos.(θ)))
    μ = period * μθ / (2π)
    return mod(μ, period)
end

function circular_diff(b, a; period=4)
    return mod(b - a + period/2, period) - period/2
end

test_offset = rx_phi_offset(t=0)

rx_phibar = mean(test_offset,dims=:ens)

for p in test_offset.path
    rx_phibar(path=p).= circular_mean(test_offset(path=p))
end

Xrx_phi = similar(test_offset)
for e in test_offset.ens
    Xrx_phi(ens=e) .= dropdims(circular_diff.(test_offset(ens=e),rx_phibar),dims=:ens)
end
