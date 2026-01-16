using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace MathTools.Core.Extensions
{
    public static class MathNetExtensions
    {
        public static Vector<double> ToVector(this double[] values)
        {
            return Vector<double>.Build.DenseOfArray(values);
        }

        public static Vector<double> ToVector(this int[] values)
        {
            var data = new double[values.Length];
            for (int i = 0; i < values.Length; i++)
            {
                data[i] = values[i];
            }

            return Vector<double>.Build.DenseOfArray(data);
        }

        public static Vector<float> ToVector(this float[] values)
        {
            return Vector<float>.Build.DenseOfArray(values);
        }

        public static double[] ToArray(this Vector<double> vector)
        {
            return vector.ToArray();
        }

        public static float[] ToArray(this Vector<float> vector)
        {
            return vector.ToArray();
        }

        public static int[] ToIntArray(this Vector<double> vector)
        {
            var data = vector.ToArray();
            var result = new int[data.Length];
            for (int i = 0; i < data.Length; i++)
            {
                result[i] = (int)data[i];
            }

            return result;
        }

        public static Matrix<double> ToMatrix(this double[][] values)
        {
            return Matrix<double>.Build.DenseOfRowArrays(values);
        }

        public static Matrix<double> ToMatrix(this int[][] values)
        {
            int rows = values.Length;
            int cols = rows == 0 ? 0 : values[0].Length;
            var data = new double[rows, cols];
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    data[i, j] = values[i][j];
                }
            }

            return Matrix<double>.Build.DenseOfArray(data);
        }

        public static Matrix<float> ToMatrix(this float[][] values)
        {
            return Matrix<float>.Build.DenseOfRowArrays(values);
        }

        public static SparseMatrix ToSparseMatrix(this double[][] values)
        {
            int rows = values.Length;
            int cols = rows == 0 ? 0 : values[0].Length;
            var entries = new System.Collections.Generic.List<Tuple<int, int, double>>();

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    double value = values[i][j];
                    if (value != 0.0)
                    {
                        entries.Add(Tuple.Create(i, j, value));
                    }
                }
            }

            return SparseMatrix.OfIndexed(rows, cols, entries);
        }

        public static SparseMatrix ToSparseMatrix(this int[][] values)
        {
            int rows = values.Length;
            int cols = rows == 0 ? 0 : values[0].Length;
            var entries = new System.Collections.Generic.List<Tuple<int, int, double>>();

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    int value = values[i][j];
                    if (value != 0)
                    {
                        entries.Add(Tuple.Create(i, j, (double)value));
                    }
                }
            }

            return SparseMatrix.OfIndexed(rows, cols, entries);
        }

        public static SparseMatrix ToSparseMatrix(this float[][] values)
        {
            int rows = values.Length;
            int cols = rows == 0 ? 0 : values[0].Length;
            var entries = new System.Collections.Generic.List<Tuple<int, int, double>>();

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    float value = values[i][j];
                    if (value != 0.0f)
                    {
                        entries.Add(Tuple.Create(i, j, (double)value));
                    }
                }
            }

            return SparseMatrix.OfIndexed(rows, cols, entries);
        }

        public static double[][] ToJaggedArray(this Matrix<double> matrix)
        {
            int rows = matrix.RowCount;
            int cols = matrix.ColumnCount;
            var result = new double[rows][];
            for (int i = 0; i < rows; i++)
            {
                result[i] = new double[cols];
                for (int j = 0; j < cols; j++)
                {
                    result[i][j] = matrix[i, j];
                }
            }

            return result;
        }

        public static float[][] ToJaggedArray(this Matrix<float> matrix)
        {
            int rows = matrix.RowCount;
            int cols = matrix.ColumnCount;
            var result = new float[rows][];
            for (int i = 0; i < rows; i++)
            {
                result[i] = new float[cols];
                for (int j = 0; j < cols; j++)
                {
                    result[i][j] = matrix[i, j];
                }
            }

            return result;
        }

        public static int[][] ToIntJaggedArray(this Matrix<double> matrix)
        {
            int rows = matrix.RowCount;
            int cols = matrix.ColumnCount;
            var result = new int[rows][];
            for (int i = 0; i < rows; i++)
            {
                result[i] = new int[cols];
                for (int j = 0; j < cols; j++)
                {
                    result[i][j] = (int)matrix[i, j];
                }
            }

            return result;
        }
    }
}
