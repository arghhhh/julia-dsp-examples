

include( "env.jl" )

import FixPts, FixPts.FixPt
import FixPts.uFixPt, FixPts.sFixPt, FixPts.split

using Test

# https://en.wikipedia.org/wiki/Karatsuba_algorithm
# The algorithm has been modified to accomodate different widths for x and y
# The FixPt types carry information on the significance of the bits in the numbers
# and because of this, the meaning of the variables below is different from the wikipedia 
# description.

# Karatsuba multiplication has a place in hardware design - but this example is not intended 
# to illustrate a useful application.
# For a more useful application, see for example:
#
# https://news.ycombinator.com/item?id=43372227
# https://arxiv.org/abs/2501.08889


# single application (not recursive)

function karatsuba_multiply( x::FixPts.FixPt, y::FixPts.FixPt, Bx = div( FixPts.bitwidth(x) + FixPts.bitwidth(y), 4 ), By = Bx )

        x1,x0 = split(x,Bx)  # x == x1 + x0
        y1,y0 = split(y,By)  # y == y1 + y0

        z2 = x1 * y1   # multiply #1
        z0 = x0 * y0   # multiply #2

        z3 = (x1 + x0<<By )*( y1>>By + y0 ) # multiply #3
        z1 = z3 - (z2 >>  By ) - (z0 << By)

        z = z2 + z1 + z0

        # note the result type has wider bounds than necessary
        # it would be possible to convert to 
        # typeof( x*y )
        return z
end

# Karatsuba multiplication is not a great example for a FixPt library
# - it is unusual to want to add different parts of numbers to each other in this way
# But it is a good example, in that things just work - for signed, unsigned, different split points
# - even split points outside of the range of valid bits (which is not a useful use-case!)

# This is also an example of a type of misuse - this code is quite slow, 
# because Julia is compiling code for every iteration for different fixed point formats
# The FixPt library is intended for situations where the variables are truely fixed point 
# - ie not changing wordlengths and binary point position.

t = []

@time @testset begin 
        for ii in 1:100
                i = rand(-10:10 )
                q = rand( -i:10 )
                x = rand( rand(Bool) ? uFixPt(i,q) : sFixPt(i,q) )
                Bx = rand(-q-1:i+1)

                i = rand(-10:10 )
                q = rand( -i:10 )
                y = rand( rand(Bool) ? uFixPt(i,q) : sFixPt(i,q) )
                By = rand(-q-1:i+1)

                push!( t, (x,y,Bx,By) )

                z1 = x * y

                @test z1 == karatsuba_multiply( x, y )
                @test z1 == karatsuba_multiply( x, y, Bx, By )
        end
end # testset

@time @testset begin 
        for ii in 1:100000
                i = rand( [-10,-5, 0, 5] )
                q = 15
                x = rand( rand(Bool) ? uFixPt(i,q) : sFixPt(i,q) )

                i = rand( [-10,-5, 0, 5 ] )
                q = 15
                y = rand( rand(Bool) ? uFixPt(i,q) : sFixPt(i,q) )

                push!( t, (x,y) )

                z1 = x * y

                @test z1 == karatsuba_multiply( x, y )
        end
end # testset
