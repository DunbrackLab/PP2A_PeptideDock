using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using PdbxmmCIFLib;

namespace PepStructInfoLib  
{
    /// <summary>
    /// Summary description for AtomInfo.
    /// </summary>
    public class AtomInfo : ICloneable
    {
        public AtomInfo()
        {
        }

        #region member variables
        public string residue = "";
        public string seqId = "";

        public int atomId = 0;
        public Coordinate xyz = new Coordinate();
        public string atomName = "";
        public string atomType = "";

        public string authSeqId = "";
        public string authResidue = "";
        public double bFac = 0.0; // added this value
        public string altConfID = ""; // added this value
        public double occupancy = 0.0; // added this value
        #endregion

        /// <summary>
        /// deep copy of external atom information
        /// </summary>
        /// <param name="extAtomInfo"></param>
        public AtomInfo(AtomInfo extAtomInfo)
        {
            this.residue = extAtomInfo.residue;
            this.seqId = extAtomInfo.seqId;
            this.authSeqId = extAtomInfo.authSeqId;
            this.atomId = extAtomInfo.atomId;
            this.authResidue = extAtomInfo.authResidue;
            this.xyz = new Coordinate(extAtomInfo.xyz.Y, extAtomInfo.xyz.Y, extAtomInfo.xyz.Z);
            this.atomName = extAtomInfo.atomName;
            this.atomType = extAtomInfo.atomType;
            this.bFac = extAtomInfo.bFac;
            this.altConfID = extAtomInfo.altConfID;
            this.occupancy = extAtomInfo.occupancy;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="cifAtom"></param>
        public AtomInfo (PdbxmmCIFAtom cifAtom)
        {
            this.residue = cifAtom.atomLineFields[5];
            this.seqId = cifAtom.atomLineFields[8];
            this.authSeqId = cifAtom.atomLineFields[16];
            this.atomId = Convert.ToInt32 (cifAtom.atomLineFields[1]);
            this.authResidue = cifAtom.atomLineFields[17];
            double x = Convert.ToDouble(cifAtom.atomLineFields[10]);
            double y = Convert.ToDouble (cifAtom.atomLineFields[11]);
            double z = Convert.ToDouble (cifAtom.atomLineFields[12]);
            this.xyz = new Coordinate(x, y, z);
            this.atomName = cifAtom.atomLineFields[3];
            this.atomType = cifAtom.atomLineFields[2];
            this.bFac = Convert.ToDouble (cifAtom.atomLineFields[14]);
            this.altConfID = cifAtom.atomLineFields[4];
            this.occupancy = Convert.ToDouble(cifAtom.atomLineFields[13]);
        }

        /// <summary>
        /// deep copy every fields
        /// </summary>
        /// <returns></returns>
        public object Clone()
        {
            // showdow copy
            AtomInfo clonedAtomInfo = (AtomInfo)this.MemberwiseClone();
            // deep copy fields with reference type
            clonedAtomInfo.xyz = new Coordinate(xyz.X, xyz.Y, xyz.Z);

            return clonedAtomInfo;
        }

        /// <summary>
        /// Euclidean distance between two 3d points
        /// </summary>
        /// <param name="atomA"></param>
        /// <param name="atomB"></param>
        /// <returns></returns>
        public static double operator -(AtomInfo atomA, AtomInfo atomB)
        {
            double sqrSum = 0.0;
            sqrSum += (atomA.xyz.X - atomB.xyz.X) * (atomA.xyz.X - atomB.xyz.X);
            sqrSum += (atomA.xyz.Y - atomB.xyz.Y) * (atomA.xyz.Y - atomB.xyz.Y);
            sqrSum += (atomA.xyz.Z - atomB.xyz.Z) * (atomA.xyz.Z - atomB.xyz.Z);
            return Math.Sqrt(sqrSum);
        }
    }
}
