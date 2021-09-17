using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;
using DbscanImplementation;
using AuxFuncLib;

namespace StructClusteringLib
{
    public class DbScanCluster
    {
        public delegate double CalculateDiAngleDist(PepDiAngleDatasetItem item1, PepDiAngleDatasetItem item2);
        public delegate double CalculateRmsdDist(PepCoordinateDatasetItem item1, PepCoordinateDatasetItem item2);

        public StructDistance structDist = new StructDistance();

        /// <summary>
        /// 
        /// </summary>
        /// <param name="diAnglesFile"></param>
        /// <param name="resultFile"></param>
        public void ClusterPeptideStructuresOnDiAngles (string diAnglesFile, string resultFile)
        {
            PepDiAngleDatasetItem[] diAngleItems = ReadDiAnglesDataFromFile(diAnglesFile, true);
            double epsilon = 0.5;
            int minPts = 20;

            HashSet<PepDiAngleDatasetItem[]> clusters = ClusterPeptideStructuresOnDiAngles(diAngleItems, epsilon, minPts);
            StreamWriter dataWriter = new StreamWriter(resultFile);
            dataWriter.WriteLine("ClusterID\tStructNo\tStructName");
            int clusterId = 1;
            foreach (PepDiAngleDatasetItem[] cluster in clusters)
            {
                foreach (PepDiAngleDatasetItem item in cluster)
                {
                    dataWriter.WriteLine(clusterId + "\t" + item.structIndex + "\t" + item.structName + "\t" +  FormatDiAnglesList (item.diAnglesList, '\t'));
                }
                clusterId++;
            }
            dataWriter.Close();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="diAnglesItems"></param>
        /// <param name="epsilon"></param>
        /// <param name="minPts"></param>
        public HashSet<PepDiAngleDatasetItem[]> ClusterPeptideStructuresOnDiAngles(PepDiAngleDatasetItem[] diAnglesItems, double epsilon, int minPts)
        {
            DbscanAlgorithm<PepDiAngleDatasetItem> dbScan = new DbscanAlgorithm<PepDiAngleDatasetItem>(CalculateDiAngleSumDist);
            HashSet<PepDiAngleDatasetItem[]> clusters = null;
            dbScan.ComputeClusterDbscan(diAnglesItems, epsilon, minPts, out clusters);

            return clusters;
        }
 
        /// <summary>
        /// 
        /// </summary>
        /// <param name="diAngleFile"></param>
        /// <param name="distMatrixFile"></param>
        /// <param name="method"></param>
        /// <returns></returns>
        public double[,] CalculateDistanceMatrix (string diAngleFile, string distMatrixFile, string method)
        {
            method = method.ToLower();
            PepDiAngleDatasetItem[] diAnglesItems = ReadDiAnglesDataFromFile(diAngleFile, true);
            double[,] distMatrix = new  double[diAnglesItems.Length, diAnglesItems.Length];
            
            for (int i = 0; i < diAnglesItems.Length; i ++)
            {
                for (int j = i + 1; j < diAnglesItems.Length; j ++)
                {
                    if (method == "sum")
                    {
                        distMatrix[i, j] = CalculateDiAngleSumDist(diAnglesItems[i], diAnglesItems[j]);
                    }
                    else if (method == "max")
                    {
                        distMatrix[i, j] = CalculateDiAngleMaxDist(diAnglesItems[i], diAnglesItems[j]);
                    }
                    distMatrix[j, i] = distMatrix[i, j];
                }
            }
            StreamWriter distWriter = new StreamWriter(distMatrixFile);
            for (int i = 0; i < diAnglesItems.Length; i++)
            {
                double[] rowItems = GetRow(distMatrix, i);
                distWriter.WriteLine(ParseHelper.FormatArrayString(rowItems, '\t'));
            }
            distWriter.Close();
            return distMatrix;
        }
    

        #region implementation for delegates
        /// <summary>
        /// 
        /// </summary>
        /// <param name="item1"></param>
        /// <param name="item2"></param>
        /// <returns></returns>
        public double CalculateDiAngleSumDist (PepDiAngleDatasetItem item1, PepDiAngleDatasetItem item2)
        {
            double sumDist = structDist.CalculateSumDihedrealAngleDist(item1, item2);
            return sumDist;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="item1"></param>
        /// <param name="item2"></param>
        /// <returns></returns>
        public double CalculateDiAngleMaxDist(PepDiAngleDatasetItem item1, PepDiAngleDatasetItem item2)
        {
            double maxDist = structDist.CalculateMaxDihedrealAngleDist(item1, item2);
            return maxDist;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="item1"></param>
        /// <param name="item2"></param>
        /// <returns></returns>
        public double CalculateAvgRmsdDist (PepCoordinateDatasetItem item1, PepCoordinateDatasetItem item2)
        {
            double avgRmsdDist = structDist.CalculateRmsdAvg(item1.coordinatesList.ToArray(), item2.coordinatesList.ToArray());
            return avgRmsdDist;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="item1"></param>
        /// <param name="item2"></param>
        /// <returns></returns>
        public double CalculateSumRmsdDist(PepCoordinateDatasetItem item1, PepCoordinateDatasetItem item2)
        {
            double sumRmsdDist = structDist.CalculateRmsd(item1.coordinatesList.ToArray(), item2.coordinatesList.ToArray());
            return sumRmsdDist;
        }
        #endregion

        #region I/O functions
        /// <summary>
        /// 
        /// </summary>
        /// <param name="clusterFiles"></param>
        /// <param name="diAnglesStructFile"></param>
        /// <returns></returns>
        public string[] AddListClustersToDataFile (string[] clusterFiles, string diAnglesStructFile)
        {
            string[] clusterAngleFiles =  new string[clusterFiles.Length];
            int fileIndex = 0;
            string clusterAddFile = "";
            foreach (string clusterFile in clusterFiles)
            {
                clusterAddFile = AddClustersToDataFile(clusterFile, diAnglesStructFile);
                clusterAngleFiles[fileIndex] = clusterAddFile;
                fileIndex++;
            }
            return clusterAngleFiles;
        }
        /// <summary>
        /// 
        /// </summary>
        public string AddClustersToDataFile (string rClusterFile, string diAnglesStructFile)
        {
            FileInfo fileInfo = new FileInfo  (diAnglesStructFile);            
            string epsilon = "";
            string minPts = "";
            string[] clusterIds = ReadClusterIds(rClusterFile, out epsilon, out minPts);

            string clusterDiAnglesFile = Path.Combine(fileInfo.DirectoryName, fileInfo.Name.Replace (".txt", "clusters_eps" + epsilon + "_minPts" + minPts + ".txt"));
            StreamWriter dataWriter = new StreamWriter(clusterDiAnglesFile);
            StreamReader dataReader = new StreamReader(diAnglesStructFile);
            string line = "";
            string headerLine = dataReader.ReadLine();
            dataWriter.WriteLine("ClusterID\t" + headerLine);
            int structIndex = 0;
            while ((line = dataReader.ReadLine ()) != null)
            {
                if (structIndex >= clusterIds.Length)
                {
                    throw new Exception("the number of cluster IDs are not equal to the number of structures. This is not correct. ");
                }
                dataWriter.WriteLine(clusterIds[structIndex] + "\t" + line);
                structIndex++;
            }
            dataReader.Close();
            dataWriter.Close();

            return clusterDiAnglesFile;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="rClusterFile"></param>
        /// <param name="epsilon"></param>
        /// <param name="minPts"></param>
        /// <returns></returns>
        public string[] ReadClusterIds (string rClusterFile, out string epsilon, out string minPts)
        {
            List<string> clusterIdList = new List<string>();
            StreamReader dataReader = new StreamReader(rClusterFile);
            string line = "";
            bool clusterStart = false;
            epsilon = "";
            minPts = "";
            int epsIndex = -1;
            int minPtsIndex = -1;
            while ((line = dataReader.ReadLine ()) != null)
            {
                if (line == "")
                {
                    continue;
                }
                epsIndex = line.IndexOf("> epsilon = ");
                minPtsIndex = line.IndexOf("> minPts = ");
                if (epsIndex > -1)
                {
                    // > epsilon = 10;
                    string[] fields = line.Split("=;".ToCharArray ());
                    epsilon = fields[1].Trim ();
                    continue;
                }
                if (minPtsIndex > -1)
                {
                    // > minPts = 5;
                    string[] fields = line.Split("=;".ToCharArray());
                    minPts = fields[1].Trim ();
                    continue;
                }
                if (line.IndexOf ("$cluster") > -1)
                {
                    clusterStart = true;
                    continue;
                }
                if (clusterStart)
                {
                    //  [1]  1  1  1  1  1  1  2  2  2  2  2  2  2  1  1  1  1  1  1  2  2
                    //  [54]  2  2  2  2  3  3  3  3  3  3  1  0  1
                    string[] fields = ParseHelper.SplitPlus(line, ' ');
                    string[] lineClusters = new string[fields.Length - 1];
                    Array.Copy(fields, 1, lineClusters, 0, lineClusters.Length);
                    clusterIdList.AddRange(lineClusters);
                }
            }
            dataReader.Close();
            return clusterIdList.ToArray();
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="diAngleFile"></param>
        /// <returns></returns>
        public PepDiAngleDatasetItem[] ReadDiAnglesDataFromFile(string diAngleFile, bool isHeader)
        {
            List<PepDiAngleDatasetItem> diAngleItemList = new List<PepDiAngleDatasetItem>();
            StreamReader dataReader = new StreamReader(diAngleFile);
            string line = "";
            string headerLine = "";
            if (isHeader)
            {
                headerLine = dataReader.ReadLine();
            }
            int structNo = 0;
            string structName = "";
            while ((line = dataReader.ReadLine()) != null)
            {
                string[] fields = line.Split('\t');
                List<double[]> diAngleList = new List<double[]>();
                structNo = Convert.ToInt32(fields[0]);
                structName = fields[1];
                for (int i = 5; i < fields.Length; i += 5)
                {
                    double[] diAngles = new double[3];
                    diAngles[0] = Convert.ToDouble(fields[i]);
                    diAngles[1] = Convert.ToDouble(fields[i + 1]);
                    diAngles[2] = Convert.ToDouble(fields[i + 2]);
                    diAngleList.Add(diAngles);
                }
                PepDiAngleDatasetItem item = new PepDiAngleDatasetItem(diAngleList, structName, structNo);
                diAngleItemList.Add(item);
            }
            dataReader.Close();
            return diAngleItemList.ToArray();
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="diAnglesList"></param>
        /// <param name="delimitor"></param>
        /// <returns></returns>
        private string FormatDiAnglesList(List<double[]> diAnglesList, char delimitor)
        {
            string diAnglesString = "";
            foreach (double[] diAngles in diAnglesList)
            {
                diAnglesString += (ParseHelper.FormatArrayString<double>(diAngles, delimitor) + delimitor);
            }
            return diAnglesString.TrimEnd(delimitor);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="matrix"></param>
        /// <param name="columnNumber"></param>
        /// <returns></returns>
        public T[] GetColumn<T> (T[,] matrix, int columnNumber)
        {
            return Enumerable.Range(0, matrix.GetLength(0))
                .Select(x => matrix[x, columnNumber])
                .ToArray();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="matrix"></param>
        /// <param name="rowNumber"></param>
        /// <returns></returns>
        public T[] GetRow<T> (T[,] matrix, int rowNumber)
        {
            return Enumerable.Range(0, matrix.GetLength(1))
                    .Select(x => matrix[rowNumber, x])
                    .ToArray();
        }
        #endregion
    }
}
