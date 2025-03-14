

include( "env.jl" )


import FixPts: FixPts, FixPt, sFixPt, uFixPt, mFixPt
import Bints: Bints, Bint
import Mints: Mints, Mint

import Processors
import Sequences
import FixedWidths

import Spectrum
import SNR

import ProcSeqs


include( "dsp-components.jl" )

x = mFixPt(5)

int1 = Processors.IIR_poles( [ 1, -1 ] )
diff1 = Processors.FIR( [ 1, -1 ] )

Sequences.concatenate( [0, 1 ], Sequences.sequence( Mint{5}(0) ), Mint{5}(0) )  |> Processors.Take(10) |> collect

Sequences.concatenate( [0, 1 ], Sequences.sequence( Mint{5}(0) ), Mint{5}(0) ) |> int1 |> Processors.Take(10) |> collect

Sequences.concatenate( [0, 1 ], Sequences.sequence( mFixPt(5,0)(0) ), mFixPt(5,0)(0) )  |> Processors.Take(10) |> collect
Sequences.concatenate( [0, 1 ], Sequences.sequence( mFixPt(5,0)(0) ), mFixPt(5,0)(0) ) |> int1 |> Processors.Take(10) |> collect

#Sequences.Sinusoid(0.01) |> Processors.MapT{ sFixPt(1,15)}( x->mFixPt{32,-16}(round( sFixPt(1,15), x )) ) 
#Sequences.Sinusoid(0.01) |> Processors.MapT{ mFixPt(32,-16) }( x->mFixPt(32,-16)(round( sFixPt(1,15), x )) ) |> Processors.Take(101) |> collect


xs = Sequences.Sinusoid(0.01) |> Processors.MapT{ mFixPt(32,-16)}( x->mFixPt(32,-16)(round( sFixPt(1,15), x )) ) 
# xs = Sequences.Sinusoid(0.01) |> Processors.MapT{ sFixPt(1,15) }( x->round( sFixPt(1,15), x ))  

xs |> Truncate_lsbs_to{-8}() |> Processors.Take(100) |> collect


xs |> int1 |> Processors.Downsample(2) |> diff1 |> Processors.Convert{ Float64 }() |> Processors.Take(100) |> collect


adc = Round{sFixPt(1,15)}()


xs = Sequences.sequence(0.01)  # DC

xs = Sequences.Sinusoid(0.01)

xs1 = xs |> Round{sFixPt(1,15)}() |> Processors.Convert{ mFixPt(32,-16) }() 

ys = xs1 |> int1 |> Processors.Downsample(20) |> Truncate_lsbs_to{0}() |> diff1 |> Processors.Convert{ Float64 }() |> Processors.Take(100) |> collect


ys = xs1 |> int1^4 |> Truncate_lsbs_to{0}() |> diff1^4 |> Processors.Convert{ Float64 }() |> Processors.Take(100000) |> collect

y1 = ys[ end-65535 : end ]

mag, residual, bhat = SNR.determine_snr( y1, [0.01] )

s = Spectrum.spectrum( residual )

using Plots

plot( cumsum(s) .|> x->10log10(x), xaxis=:log, legend=nothing )
plot!( [10,100,1000,10000],[-90*3,-90*2,-90,0] )


@show sum(s)

ntf = diff1^4
f = range( 0.0001, 0.5 ; length = 1001 )

h = ProcSeqs.freqz2( ntf, f )

bw = 1 / length(f)  # approx

cummulative_noise = cumsum( h ) / length(f) / 12
@show  sum( h ) / length(f) / 12

plot!( f .* 65536, h .|> x->10log10(x) )
plot!( f .* 65536, cummulative_noise .|> x->10log10(x) )

ys = rand( Float64, 65536 ) .- 0.5
s = Spectrum.spectrum( ys )
sum(s)

f = range( 0.0001, 0.5 ; length = 1001 )
h = ones( length(f) )
cummulative_noise = cumsum( h ) / length(f) / 12



ntf = diff1
h = ProcSeqs.freqz2( ntf, f )
@show cummulative_noise = sum( h ) / length(f) / 12

rands = ys = rand( Float64, 65536 ) .- 0.5

ys = rands |> ntf |> collect
s = Spectrum.spectrum( ys )
@show sum(s)






