
# Constantinides, Cheung & Luk paper (Electronics Letters 11th Nov 1999, Vol 35, No.23)
# ----------------------------------

# Apparently the LSB size / 12 quantization noise "rule" is not correct when 
# applied to requantization within digital signal processing systems...
#
# The assumptions leading to 1/12 rule include the noise being continuously distributed 
# over error range - whereas in the requantization situation, the errors are discrete 
# valued.  The mean and variance for the discrete case are easily derived using the 
# the sum of an arithmetic series and the sum of squares of an arithmetic series.
#
# This is nothing to do with correlation between the input and error - which also 
# needs to be considered.
#

# not using julia-fixed-point - not needed here
# not using julia-signal-systems - very simple and memoryless system


using Plots

import Statistics

quant(x) = floor(x)
quant_error(x) = quant(x) - x

scale = 1000000

num_samples = 100000000

xs = scale * rand( num_samples )


xs = range( -5, 5 ; length = 1000 )
ys = quant_error.(xs)

plot( xs, ys )

# expect this to be -0.5 (output is rounded down)
@show Statistics.mean( ys )
# expect this to be 1/12  (0.0833333333)
@show Statistics.var( ys )


# Now look at situation where the input is already quantized
# Consider the case of the input having N quantization levels

# Using truncation of an unsigned or twos-complement number
# the error is unformly drawn from -(0..N-1)

N = 8
errors = -(0:N-1)/N

ys = rand( errors, num_samples )

# expecting -1/2 * (N-1)/N
@show -1/2 * (N-1)/N
@show Statistics.mean( ys )

# expecting 1/12*(N^2-1)/(N^2)
@show 1/12*(N^2-1)/(N^2)
@show Statistics.var( ys )


requant_error_mean(N) =  -1/2 * (N-1)/N
requant_error_var(N)  =  1/12*(N^2-1)/(N^2)

Ns = 1:64

plot( log2.(Ns), [ requant_error_var.(Ns) requant_error_mean.(Ns)  ], label=nothing , seriescolor=[ :blue :red ])

Ns1 = [1,2,4,8,16,32,64]
scatter!( log2.(Ns1), [ requant_error_var.(Ns1) requant_error_mean.(Ns1)  ],  seriescolor=[ :blue :red ], label = [  "var" "mean" ] )
plot!( legend=:right )


# Oppenheim & Schafer 1975 book, p411 Eqn (9.3)
# Also p412 top 
#   2^-b1 << 2^-b  so 2^-b1 can be ignored

# Constantinides, Cheung & Luk paper (Eletcronics Letters 11th Nov 1999, Vol 35, No.23)
# They take the bounds for the truncation error developed in the Oppenheim & Schafer book
# and applies it to the uniformly continuously distributed variance expression (1/12*lsb)
# also developed in that book - but I could not see anywhere where Oppenheim & Schafer 
# did this (which would have been error)
requant_var_unsimplified(N) = 1/12*(1-1/N)
requant_var_simplified(N)   = 1/12

plot!( log2.(Ns), [ requant_var_unsimplified.(Ns) requant_var_simplified.(Ns)  ], label=nothing )
plot!( title="Error Mean & Variance after requantization" )
plot!( xlabel="Number of bits truncated" )
plot!( ylabel="Mean & variance relative to the output LSB size" )

savefig( "requantization.svg" )
