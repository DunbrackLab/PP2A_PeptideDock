using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PeptideDockingLib
{
    public class MotifSeqInfo
    {
        public string motif;   // sequence
        public string chainId;
        public int motifBeg;  // this is the sequential residue number 1...N
    }
}
