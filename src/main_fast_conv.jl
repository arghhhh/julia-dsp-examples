include( "env.jl" )

import Processors
import Sequences
import FFTW

import DSPfns

function fast_convolution( x, h )
        # x already assumed to be sufficiently padded

        h1 = zeros( eltype(h), length( x ) )
        h1[1:length(h)] .= h
        H = FFTW.rfft(h1)
        X = FFTW.rfft(x)
        Y = H .* X
        y = FFTW.irfft( Y, length(x) )
        return y
end

function fast_convolution_part1( x, h )
        # x already assumed to be sufficiently padded

        h1 = zeros( eltype(h), length( x ) )
        h1[1:length(h)] .= h
        H = FFTW.rfft(h1)

        return H
end
function fast_convolution_part2( x, H )

        X = FFTW.rfft(x)
        Y = H .* X
        y = FFTW.irfft( Y, length(x) )
        return y
end


struct FastConvolution <: Processors.SampleProcessor
        h
end
# Processors.process( p::FastConvolution, x, state=nothing ) = begin
#         fast_convolution(x,p.h),nothing
# end

Processors.process( p::FastConvolution, x ) = begin
        state = fast_convolution_part1( x, p.h )
        Processors.process( p, x, state )
end
Processors.process( p::FastConvolution, x, state ) = begin
        fast_convolution_part2(x,state),state
end

Base.eltype( ::Type{ Processors.Apply{I,FastConvolution} }) where {I} = Vector{ Float64 }





x = [ 1,2,3,4,5,6,7,8,9,10,11,12 ]
h = [1,2,3,4]

@show fast_convolution( x, h )
@show  DSPfns.conv(x,h )

n = 4


sys = (
           Processors.MapT{ Vector{Float64} }( x->fill(x,n) ) 
        |> Processors.Delays1()
        |> Processors.Downsample(n,n-1)
        |> Processors.Upsample(n)
        |> Processors.Delays2()
        |> Processors.MapT{Float64}( sum )
        )

@show 1:20 |> sys |> collect


# This:
#    (  Processors.MapT{ Vector{Int} }( x->fill(x,L+M-1) ) 
#    |> Processors.Delays1()
#    )
# is equivalent to Processors.Vectorize(L+M-1)
x = 1:100
n = 4
sys1 = Processors.MapT{ Vector{Int} }( x->fill(x,n) ) |> Processors.Delays1()
sys2 = Processors.Vectorize(n)

y1 = x |> sys1 |> collect
y2 = x |> sys2 |> (x->collect(Vector,x))

@assert y1 == y2

# @error("stop")


x = x = 1:100 |> Processors.MapT{ Vector{Float64} }( x->fill(x,n) ) |> Processors.Delays1() |> Processors.Downsample(n,n-1)

L = n

sys1 = Processors.Upsample(L) |> Processors.Delays2() |> Processors.MapT{Float64}( sum )
sys2 = Processors.Serialize()

y1 = x |> sys1 |> collect
y2 = x |> sys2 |> collect

@assert y1 == y2




# overlap-save:

# sliding window min length L+M-1 
# down sample by L
# convolve with h length M, to get result length L+2M-2 which can be considered
#   to have M-1 "bad" samples at either end, and L good samples in the middle
# reconstruct using the L good samples

h = [1,2,3,4]
h = rand( 1:10, 50 )

L = 8
M = length(h)

impulse = [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

conv_core( x, h ) = begin
        y1 = conv(x,h)
        y = y1[M:M+L-1]
        println()
        @show x h y1 y
        println()
        return y
end


overlapsave( vlen, decimation ) = Processors.Vectorize(vlen) |> Processors.Downsample(decimation,decimation-1)


sys1 = (
           overlapsave( L+M-1, L )
        |> Processors.MapT{Vector{Int}}( x->DSPfns.conv(x,h) )
        |> Processors.MapT{Vector{Int}}( x->x[M:M+L-1] )
        |> Processors.Serialize()
        )
sys2 = (
        overlapsave( L+M-1, L )
        |> FastConvolution(h)
        |> Processors.MapT{Vector{Int}}( x->x[M:M+L-1] )
        |> Processors.Serialize()
        )
xs = rand( 0:10, 2000 )

y0 = DSPfns.conv( xs, h )

y1 = xs |> sys1 |> collect
y2 = xs |> sys1 |> collect

@show length(y1) M L
@assert y1 == y0[1:length(y1)]
@assert y2 == y1

# error("stop")
# 
# sys1 = (
#            Processors.MapT{ Vector{Int} }( x->fill(x,L+M-1) ) 
#         |> Processors.Delays1()
#      #   |> Processors.Downsample(L,L-1)
# )
# 
# xs = 1:20
# xs |> sys1 |> collect


# overlap-add

# SlidingWindow(L)
# Downsample(L)
# convolve
# v->reshape(v,L)
# Delays2
# v->sum(v)
# serialize 

overlapadd( vlen ) = Processors.Reshape(vlen) |> Processors.Delays2() |> Processors.Sum()

zeropad( v, n ) = begin
        y = zeros( eltype(v), length(v) + n )
        y[1:length(v)] .= v
        return y
end

sys1 = 
        (  Processors.SlidingWindow(L) 
        |> Processors.Downsample(L) 
        |> Processors.MapT{Vector{Int}}( x->DSPfns.conv(x,h) )
    #    |> Reshape(L)
    #    |> Processors.Delays2()
    #    |> Sum()
        |> overlapadd(L)
        |> Processors.Serialize()
        )
sys2 = 
        (  Processors.SlidingWindow(L) 
        |> Processors.Downsample(L) 
        |> Processors.MapT{Vector{Int}}( x->zeropad(x,2M-2) )
        |> FastConvolution(h)
    #    |> Reshape(L)
    #    |> Processors.Delays2()
    #    |> Sum()
        |> overlapadd(L)
        |> Processors.Serialize()
        )
y0 = DSPfns.conv( xs, h )

y1 = xs |> sys1 |> collect
y2 = xs |> sys2 |> collect

@show length(y1) M L
@assert y1 == y0[1:length(y1)]
#@assert y2 == y1
@assert maximum( abs.( y1 - y2 ) ) < 1e-10
        


    #    xs |> sys |> Sequences.info



nothing


