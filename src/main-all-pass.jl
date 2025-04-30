include( "env.jl")

using Test
using Plots
import Statistics

import Sequences
import Processors
import ProcSeqs
import SNR
import Spectrum

import DSPfns
import FixPts: FixPts, FixPt, sFixPt, truncate_lsbs_to
import Processors: SampleProcessor, IIR, Delay, Downsample, process, Apply, Take, Convert
import Sequences: Sinusoid, sequence

include( "dsp-components.jl" )


# "Digital signal processing schemes for efficient interpolation and decimation"
# Valenzula and Constantinides
# IEE Proc Vol 130, Pt G, No. 6 December 1983, p225
# Fig 15, RHS

a1 = 0.157606
a2 = 0.614840

allpass2(a) = IIR( [a, 0, 1 ] , [ 1, 0, a ] )

lowpass = ( allpass2( a1 ) + ( allpass2( a2 ) |> Delay(1) ) ) / 2

f_array = range( 0.001, 0.499; length = 500 )
lowpass_res  = ProcSeqs.freqzdB( lowpass , f_array )
plot( f_array, lowpass_res, label="Floating point coefficients", linewidth = 3 )
plot!(dpi=400)
savefig( "lowpass-response-unlabeled.png" )

for q in 4:8
        a1q = round( sFixPt(1,q), a1 )
        a2q = round( sFixPt(1,q), a2 )

        local lowpass = ( allpass2( a1q ) + ( allpass2( a2q ) |> Delay(1) ) ) / 2 

        plot!( f_array, ProcSeqs.freqzdB( lowpass , f_array ),label="$q coefficient bits" )
end

plot!( legend=:bottomleft, ylim=(-80,0) )
plot!( xlabel="Frequency (fs=1)", ylabel="Response (dB)", dpi=400 )
savefig( "lowpass-freq-response.png")

plot(  f_array, ProcSeqs.freqzdB( Downsample(2)^0 |> lowpass, f_array ), label="1" )
plot!( f_array, ProcSeqs.freqzdB( Downsample(2)^1 |> lowpass, f_array ), label="2" )
plot!( f_array, ProcSeqs.freqzdB( Downsample(2)^2 |> lowpass, f_array ), label="3" )
plot!( xlabel="Frequency (fs=1)", ylabel="Response (dB)", dpi=400 )
savefig( "lowpass-x3-freq-response-a.png")

plot(  f_array, ProcSeqs.freqzdB( ( lowpass |> Downsample(2) )^1 , f_array ), label="1" )
plot!( f_array, ProcSeqs.freqzdB( ( lowpass |> Downsample(2) )^2 , f_array ), label="2" )
plot!( f_array, ProcSeqs.freqzdB( ( lowpass |> Downsample(2) )^3 , f_array ), label="3" )
plot!( xlabel="Frequency (fs=1)", ylabel="Response (dB)", dpi=400 )
savefig( "lowpass-x3-freq-response.png")

plot(  f_array, ProcSeqs.freqzdB( lowpass, f_array ), legend=nothing, xlim=(0,0.25), ylim=(-1e-4,1e-4) )
plot!( xlabel="Frequency (fs=1)", ylabel="Response (dB)", dpi=400 )
savefig( "lowpass-freq-response-detail.png")



a1q = round( sFixPt(1,10), a1 )
a2q = round( sFixPt(1,10), a2 )
lowpassq = ( allpass2( a1q ) + ( allpass2( a2q ) |> Delay(1) ) ) / 2

# -----------------------------------------------------------------
# define an implementation of the allpass structure
# that uses floating point arithmetic - this is a stepping
# point towards getting a Fixed Point implementation
struct Allpass_Float <: SampleProcessor
        a::Float64
end
function process( p::Allpass_Float, x, state=[ 0.0, 0.0 ] )
        multout = p.a * (x + state[1])
        next_state = [ state[2], x - multout ]
        y = 0.5 * ( state[1] + multout )
        return y,next_state
end
Base.eltype( ::Type{ Apply{I,Allpass_Float} }) where {I} = Float64

# -----------------------------------------------------------------

# Now define another version that uses Fixed point arithmetic
# Word lengths are not optimized.
# For a decimator, should run the allpass filters at half the rate
# and modify this appropriately - though this should not affect the results.
# Because of the signal flow graph loops, it is necessary to truncate 
# and clamp to constrain word length growth
struct Allpass_FixPt <: SampleProcessor
        a::FixPt
end
function process( p::Allpass_FixPt, x, state=sFixPt(2,15)[ 0, 0 ] )
        multout = truncate_lsbs_to(-15, p.a * (x + state[1]) )
        next_state = [ state[2], clamp( x - multout, sFixPt(2,15) ) ]
        y = state[1] + multout

        # the shift right below is to implement the 0.5 factor required for when combining two allpass filters
        # to get the low pass filter
        return sFixPt(2,16)( y>>1 ),next_state
end
Base.eltype( ::Type{ Apply{I,Allpass_FixPt} }) where {I} = sFixPt(2,16)







step_input = sequence(1.0)



ys1 = step_input |> allpass2(       a1 ) |> Take(10) |> collect
ys2 = step_input |> Allpass_Float(  a1 ) |> Take(10) |> collect

adc = Round{sFixPt(1,15)}()

fs = 1.0

xs = Sinusoid(0.1,fs) |> adc |> Take(10)

step_input2 = step_input |> adc
ys3 = step_input2 |> Allpass_FixPt(  round( sFixPt(1,7), a1 ) ) |> Take(10) |> collect

lowpass_FixPt = Allpass_FixPt( a1q ) + ( Allpass_FixPt( a2q ) |> Delay(1) )
ys4 = step_input  |> lowpass       |> Take(10) |> collect
ys5 = step_input2 |> lowpass_FixPt |> Take(10) |> collect






n_ignore  = 1000
n_samples = 16384

for fsig in ( 0.01, 0.012345 )
        ys = Sinusoid(fsig,fs) |> adc |> lowpass_FixPt |> Convert{Float64}() |> Take( n_ignore + n_samples ) |> collect
        ys = ys[n_ignore+1:end]
        mag, residual, bhat = SNR.determine_snr( ys, [fsig] ) 
        f_axis, s_sig   = Spectrum.axis_spectrum_dB( ys       )
        f_axis, s_noise = Spectrum.axis_spectrum_dB( residual )

        # plot everything:
        plot( f_axis, [ s_sig, s_noise ], ylim=(-150,0), label=[ "signal" "residual" ] )
        plot!( xlabel="Frequency (fs=1)", ylabel="dBFs" )
        plot!( title="fsig = $fsig", dpi=400 )
        savefig( "spectrum-fsig-$(fsig).png" )

        # zoom in on detail
        plot!( ylim=(-150,0), xlim=(0,0.03) )
        savefig( "spectrum-fsig-$(fsig)-detail.png" )
end




res = []

for fsig = f_array

        ys = Sinusoid(fsig,fs) |> adc |> lowpass_FixPt |> Convert{Float64}() |> Take( n_ignore + n_samples ) |> collect

         # least squares fit, with one sine component:
        mag, residual, bhat = SNR.determine_snr( ys[n_ignore+1:end], [fsig] ) 

        s = Spectrum.spectrum(residual)
        wideband_noise_est = Spectrum.bandlimited_power_dB(s, 0.0, 0.5 , fs )
        inband_noise_est   = Spectrum.bandlimited_power_dB(s, 0.0, 0.25, fs )
        y_sigdB = 20.0 * log10( mag[1] )

        # save results onto the output array as a named tuple:
        push!( res, (; wideband_noise_est, inband_noise_est , y_sigdB ) )
end


plot(  f_array, [ r.y_sigdB for r in res ], label="signal" )
plot!( f_array, [ r.wideband_noise_est for r in res ], label = "wideband noise" )
plot!( f_array, [ r.inband_noise_est for r in res ], label = "inband noise" )
plot!( f_array, ProcSeqs.freqzdB( lowpassq, f_array ),   color=:black, label="transfer function, quantized coefficients" )
plot!( xlabel="Frequency (fs=1)", ylabel="Response (dB)", dpi=400 )
plot!( legend=:left )
savefig( "lowpass-freq-sweep.png")

plot(  f_array, [ r.y_sigdB for r in res ], linewidth = 3, label="simulated & measured" )
plot!( f_array, ProcSeqs.freqzdB( lowpass , f_array ), label="transfer function" )
plot!( f_array, ProcSeqs.freqzdB( lowpassq, f_array ),   color=:black, label="transfer function, quantized coefficients" )
plot!( xlim=(0,0.2), ylim=(-5e-5,1e-5), legend=:bottomleft )
plot!( xlabel="Frequency (fs=1)", ylabel="Response (dB)", dpi=400 )
savefig( "lowpass-freq-sweep-detail.png")

