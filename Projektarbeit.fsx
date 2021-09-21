#r "nuget: FSharpAux"
#r "nuget: FSharpAux.IO"
#r "nuget: FSharp.Stats, 0.4.1" 
#r "nuget: BioFSharp, 2.0.0-beta6"
#r "nuget: BioFSharp.IO, 2.0.0-beta6"
#r "nuget: Plotly.NET, 2.0.0-beta9"

open FSharpAux
open FSharpAux.IO
open BioFSharp
open BioFSharp.IO
open FSharp.Stats
open Plotly.NET

open BioFSharp

open System
Math.Log

Math.Log10 5.

//let yourDistance (seqA: seq<IBioItem>) (seqB:seq<IBioItem>) = 
let yourSequences 
//reconstruct a phylogenetic tree from tagged sequences
let myTree =
    PhylogeneticTree.ofTaggedBioSequences
       // yourDistance // your distance function for either p, JC69, or K81 distance
        //yourSequences // your adequate nucleotide test sequences as tagged sequences 
 
open ClustalOWrapper
let cw = ClustalOWrapper("/Users/pride_exe/Desktop/clustalo")

let sequences = 
    [
    TaggedSequence.create "seq1" ("ATGAAAATTGGCAATCCTTT")
    TaggedSequence.create "seq2" ("TTTGGUAAATGUUGATGGAA")
    TaggedSequence.create "seq3" ("GUUAAAGTTTCCATTGAAUT")
    TaggedSequence.create "seq4" ("ATGGAATTCAAATAAGGTTA")
    TaggedSequence.create "seq5" ("TGGATAGAATGGAAACCUAA")
    TaggedSequence.create "seq6" ("TGGUUUTAAAGAAAAAATUG")
    TaggedSequence.create "seq7" ("GUUGAAATTAACCCTAATTT")
    TaggedSequence.create "seq8" ("UUUGAAGTCCCTAAAAGTUU")
    TaggedSequence.create "seq9" ("CTAGCATGCATGCATGCATG")
    TaggedSequence.create "seq10" ("GATCGATCGATCGATCGATC")
    ]

//Transititon = A <-> G / T <-> C 
//Transversion = A <-> T / A <-> C / G <-> T / G <-> C

let sumup (x:float) y =
    x + y
sumup 3.1 4.5
 
let seg1 =  BioList.ofNucleotideString "ATGAAAATTGGCAATCCTTT"
let seg2 =  BioList.ofNucleotideString "TTTGGUAAATGUUGATGGAA" 
            // Herausfinden ob foor-loop bei Bio.List funktioniert

let list = [1.;2.;3.;4.;5.]
for i in list do
   printfn "%f" i // %f für float ansonsten funktioniert "i" nicht -> %_ für Nucleotides.Nucleotid finden!

for i in seg1 do
   printfn "%A" i //kein Modul gefunden.. 1. Etwas finden womit %s funktioniert oder 2. eine Funktion schreiben die Bio.List in stringList zusätzlich umschreibt

let counttransition (x: Nucleotides.Nucleotide) (y: Nucleotides.Nucleotide) = //tuple
    if "A" -> "G" then add 1 //matchen -> true soll zurückgegeben werden
    elif "G" -> "A" then add 1
    elif "T" -> "C" then add 1
    elif "C" -> "T" then add 1
    else add 0
    

let alignedSequences = 
    cw.AlignSequences(sequences,Seq.empty)


let pairwiseDistance (s:float) v l = 
    (s + v) / l 

let p = pairwiseDistance 4. 7. 9.

let jc (p:float) = 
    -(3./4.)*Math.Log (1. - 4.) /3. * p
jc p

let p2 (s:float) l = 
    s / l

let q (v:float) l =
    v / l

let kimura p2 q =
    -0.5*Math.Log10 (1.-2.*p2 - q) - 0.25*Math.Log (1-2*q2)


hello