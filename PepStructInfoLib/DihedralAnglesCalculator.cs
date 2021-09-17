using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace PepStructInfoLib
{
    /// <summary>
    /// this class is modified from Benjamin North C# class PeptideAnalysis.cs
    /// 11/18/2019
    /// </summary>
    public class DihedralAnglesCalculator
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="chain"></param>
        /// <param name="atomIndex"></param>
        /// <returns></returns>
       public double[] CalculatePhiPsiAndOmega(AtomInfo[] chain, int atomIndex)
        {
            double[] phiPsiAndOmega = new double[3] { 999, 999, 999 };

            if (chain[atomIndex].atomName != "CA") // this is the calpha atom index for the residue
            {
                return phiPsiAndOmega;
            }

            int[] atomIndices = new int[6] { -1, -1, -1, atomIndex, -1, -1 };

            //order: 
            // [0] CA (index-1)
            // [1] C' (index-1)
            // [2] N
            // [3] CA
            // [4] C'
            // [5] N (index+1)

            // now, identify all the atoms necessary

            for (int backIndex = atomIndex - 1; backIndex > -1; backIndex--)
            // this would be prettier by using residue index #s but "seqID" is a string, not an int
            {
                if (chain[backIndex].atomName == "CA" && (chain[backIndex].altConfID == "" || chain[backIndex].altConfID == "A"))
                {
                    atomIndices[0] = backIndex;
                    break; // found CA for previous atom
                }
            }

            for (int backIndex = atomIndex - 1; backIndex > -1; backIndex--)
            {
                if (chain[backIndex].atomName == "CA" && (chain[backIndex].altConfID == "" || chain[backIndex].altConfID == "A"))
                {
                    break; // this means that the previous residue was hit before atoms could be found, i.e. missing atoms
                }

                if (chain[backIndex].atomName == "C" && (chain[backIndex].altConfID == "" || chain[backIndex].altConfID == "A"))
                {
                    atomIndices[1] = backIndex;
                    break; // found the C-atom, can exit now
                }
            }

            for (int backIndex = atomIndex - 1; backIndex > -1; backIndex--)
            {
                if (chain[backIndex].atomName == "CA" && (chain[backIndex].altConfID == "" || chain[backIndex].altConfID == "A"))
                {
                    break;
                }

                if (chain[backIndex].atomName == "N" && (chain[backIndex].altConfID == "" || chain[backIndex].altConfID == "A"))
                {
                    atomIndices[2] = backIndex;
                    break; // found the N-atom, can exit now
                }
            }

            for (int forwardIndex = atomIndex + 1; forwardIndex < chain.Length; forwardIndex++)
            {
                if (chain[forwardIndex].atomName == "CA" && (chain[forwardIndex].altConfID == "" || chain[forwardIndex].altConfID == "A"))
                {
                    break;
                }

                if (chain[forwardIndex].atomName == "C" && (chain[forwardIndex].altConfID == "" || chain[forwardIndex].altConfID == "A"))
                {
                    atomIndices[4] = forwardIndex;
                    break;
                }
            }

            for (int forwardIndex = atomIndex + 1; forwardIndex < chain.Length; forwardIndex++)
            {
                if (chain[forwardIndex].atomName == "CA" && (chain[forwardIndex].altConfID == "" || chain[forwardIndex].altConfID == "A"))
                {
                    break;
                }

                if (chain[forwardIndex].atomName == "N" && (chain[forwardIndex].altConfID == "" || chain[forwardIndex].altConfID == "A"))
                {
                    atomIndices[5] = forwardIndex;
                    break;
                }
            }

            if ((atomIndices[1] > -1) && (atomIndices[3] > -1))
            {
                double[] atom1 = new double[3];
                double[] atom2 = new double[3];
                double[] atom3 = new double[3];

                atom1[0] = chain[atomIndices[1]].xyz.X;
                atom1[1] = chain[atomIndices[1]].xyz.Y;
                atom1[2] = chain[atomIndices[1]].xyz.Z;

                atom2[0] = chain[atomIndices[2]].xyz.X;
                atom2[1] = chain[atomIndices[2]].xyz.Y;
                atom2[2] = chain[atomIndices[2]].xyz.Z;

                atom3[0] = chain[atomIndices[3]].xyz.X;
                atom3[1] = chain[atomIndices[3]].xyz.Y;
                atom3[2] = chain[atomIndices[3]].xyz.Z;

                if ((atomIndices[0] > -1))
                {
                    double[] atom0 = new double[3];
                    atom0[0] = chain[atomIndices[0]].xyz.X;
                    atom0[1] = chain[atomIndices[0]].xyz.Y;
                    atom0[2] = chain[atomIndices[0]].xyz.Z;
                    phiPsiAndOmega[DihedralAngleName.omega] = CalculateTorsion(atom0, atom1, atom2, atom3);   // omega
                }
                if ((atomIndices[4] > -1) && (atomIndices[5] > -1))
                {
                    double[] atom4 = new double[3];
                    double[] atom5 = new double[3];
                    atom4[0] = chain[atomIndices[4]].xyz.X;
                    atom4[1] = chain[atomIndices[4]].xyz.Y;
                    atom4[2] = chain[atomIndices[4]].xyz.Z;
                    atom5[0] = chain[atomIndices[5]].xyz.X;
                    atom5[1] = chain[atomIndices[5]].xyz.Y;
                    atom5[2] = chain[atomIndices[5]].xyz.Z;
                    phiPsiAndOmega[DihedralAngleName.phi] = CalculateTorsion(atom1, atom2, atom3, atom4); // phi
                    phiPsiAndOmega[DihedralAngleName.psi] = CalculateTorsion(atom2, atom3, atom4, atom5);  // psi
                }
            }
            return phiPsiAndOmega;
        }

        /// <summary>
       ///  //order: 
       // [0] CA (index-1)
       // [1] C' (index-1)
       // [2] N
       // [3] CA
       // [4] C'
       // [5] N (index+1)
        /// </summary>
        /// <param name="angleAtoms">must be in the order</param>
        /// <returns>in the order of omega, phi, psi</returns>
        public double[] CalculatePhiPsiAndOmega (AtomInfo[] angleAtoms)
        {
            double omega = CalculateTorsion(angleAtoms[0].xyz, angleAtoms[1].xyz, angleAtoms[2].xyz, angleAtoms[3].xyz);
            double phi = CalculateTorsion(angleAtoms[1].xyz, angleAtoms[2].xyz, angleAtoms[3].xyz, angleAtoms[4].xyz);
            double psi = CalculateTorsion(angleAtoms[2].xyz, angleAtoms[3].xyz, angleAtoms[4].xyz, angleAtoms[5].xyz);
            double[] dihedralAngles = new double[3];
            dihedralAngles[DihedralAngleName.omega] = omega;
            dihedralAngles[DihedralAngleName.phi] = phi;
            dihedralAngles[DihedralAngleName.psi] = psi;
            return dihedralAngles;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="_atom1"></param>
        /// <param name="_atom2"></param>
        /// <param name="_atom3"></param>
        /// <param name="_atom4"></param>
        /// <returns></returns>
        public double CalculateTorsion(Coordinate atom1, Coordinate atom2, Coordinate atom3, Coordinate atom4)
        {
            double[] vec1 = new double[3] { atom2.X - atom1.X, atom2.Y - atom1.Y, atom2.Z - atom1.Z };
            double[] vec2 = new double[3] { atom3.X - atom2.X, atom3.Y - atom2.Y, atom3.Z - atom2.Z };
            double[] vec3 = new double[3] { atom4.X - atom3.X, atom4.Y - atom3.Y, atom4.Z - atom3.Z };

            double[] projVec = new double[3] { Math.Sqrt(DotProduct(vec2, vec2)) * vec1[0], Math.Sqrt(DotProduct(vec2, vec2)) * vec1[1],
                Math.Sqrt(DotProduct(vec2, vec2)) * vec1[2] };

            return (Math.Atan2((DotProduct(projVec, CrossProduct(vec2, vec3))),
                DotProduct(CrossProduct(vec1, vec2), CrossProduct(vec2, vec3)))) * 180 / Math.PI;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_atom1"></param>
        /// <param name="_atom2"></param>
        /// <param name="_atom3"></param>
        /// <param name="_atom4"></param>
        /// <returns></returns>
        public double CalculateTorsion(double[] atom1, double[] atom2, double[] atom3, double[] atom4)
        {
            double[] vec1 = new double[3] { atom2[0] - atom1[0], atom2[1] - atom1[1], atom2[2] - atom1[2] };
            double[] vec2 = new double[3] { atom3[0] - atom2[0], atom3[1] - atom2[1], atom3[2] - atom2[2] };
            double[] vec3 = new double[3] { atom4[0] - atom3[0], atom4[1] - atom3[1], atom4[2] - atom3[2] };

            double[] projVec = new double[3] { Math.Sqrt(DotProduct(vec2, vec2)) * vec1[0], Math.Sqrt(DotProduct(vec2, vec2)) * vec1[1],
                Math.Sqrt(DotProduct(vec2, vec2)) * vec1[2] };

            return (Math.Atan2((DotProduct(projVec, CrossProduct(vec2, vec3))),
                DotProduct(CrossProduct(vec1, vec2), CrossProduct(vec2, vec3)))) * 180 / Math.PI;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_vec1"></param>
        /// <param name="_vec2"></param>
        /// <returns></returns>
        private double[] CrossProduct(double[] vec1, double[] vec2)
        {
            double[] crossProductValue = new double[3];
            crossProductValue[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
            crossProductValue[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
            crossProductValue[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
            return crossProductValue;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_vec1"></param>
        /// <param name="_vec2"></param>
        /// <returns></returns>
        private double DotProduct(double[] vec1, double[] vec2)
        {
            return (vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2]);
        }
        
    }
}
