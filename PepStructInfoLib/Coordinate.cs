using System;
using System.Collections.Generic;

namespace PepStructInfoLib
{
    public class Coordinate
    {
        public double X;
        public double Y;
        public double Z;

        public Coordinate ()
        {
            X = 0;
            Y = 0;
            Z = 0;
        }

        public Coordinate (double x, double y, double z)
        {
            this.X = x;
            this.Y = y;
            this.Z = z;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="coord1"></param>
        /// <param name="coord2"></param>
        /// <returns></returns>
        public static double operator - (Coordinate coord1, Coordinate coord2)
        {
            double sqrSum = 0.0;
            sqrSum += Math.Pow (coord1.X - coord2.X, 2);
            sqrSum += Math.Pow(coord1.Y - coord2.Y, 2);
            sqrSum += Math.Pow(coord1.Z - coord2.Z, 2);

            return Math.Sqrt(sqrSum);
        }
    }
}
