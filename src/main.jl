

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

Sequences.concatenate( FixPt{0}[0, 1 ], Sequences.sequence( mFixPt(5,0)(0) ), mFixPt(5,0)(0) )  |> Processors.Take(10) |> collect
Sequences.concatenate( FixPt{0}[0, 1 ], Sequences.sequence( mFixPt(5,0)(0) ), mFixPt(5,0)(0) ) |> int1 |> Processors.Take(10) |> collect

#Sequences.Sinusoid(0.01) |> Processors.MapT{ sFixPt(1,15)}( x->mFixPt{32,-16}(round( sFixPt(1,15), x )) ) 
#Sequences.Sinusoid(0.01) |> Processors.MapT{ mFixPt(32,-16) }( x->mFixPt(32,-16)(round( sFixPt(1,15), x )) ) |> Processors.Take(101) |> collect






