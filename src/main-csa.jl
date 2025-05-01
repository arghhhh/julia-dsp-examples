

include( "env.jl" )

import FixPts, FixPts.FixPt
import FixPts.uFixPt, FixPts.sFixPt, FixPts.split

using Test


# this is a simple test case showing that the two parts of the carry save result
# can be considered as bounded integers (and not modulo 2^N)

# it may be possible to shorten these wordlengths, make it more general etc

# The main point is, the carry,save numbers can be sign extended - and so can 
# be shifted right with out error.



function csa_add( a::FixPt, b::FixPt, c::FixPt )

        Ta = typeof(a)
        Tb = typeof(b)
        Tc = typeof(c)
 
        # type level:
        Tresult = Ta + Tb + Tc   

        Tresult = promote_type( typeof(a), typeof(b), typeof(c) )

        a1 = Tresult(a)
        b1 = Tresult(b)
        c1 = Tresult(c)

        w = FixPts.bitwidth( Tresult )
        msb_mask = 1 << (w-1)

        # go from FixPt => Bint => int
        an = a1 |> Integer |> Integer
        bn = b1 |> Integer |> Integer
        cn = c1 |> Integer |> Integer

        # do the arithmetic - no carry chains
        # these are independent (parallel) bit-wise operations:
        s_n = xor( xor( an, bn ), cn )
        co_n = an & bn | an & cn | bn & cn
        
        # the carry out has the same width - but twice the significance
        # so make result FixPt types and shift the carry out left one place:
        co = Tresult( co_n ) << 1
        s  = Tresult( s_n )

        return co,s
end


ibits = 6
r = -32:31

a = sFixPt(ibits,0)(7)
b = sFixPt(ibits,0)(-8)
c = sFixPt(ibits,0)(-3)

csa_add(a,b,c)

# TODO: add iteration of Bint and FixPt types
# eg for x in sFixPt(3,2)

@testset begin
        for an in r
                a = sFixPt(ibits,0)(an)
                
                for bn in r
                        b = sFixPt(ibits,0)(bn)

                        for cn in r
                                c = sFixPt(ibits,0)(cn)

                                co,s = csa_add( a, b, c )

                                y = co + s

                                @test an+bn+cn == y.n.n
                        end
                end
        end
end
nothing
