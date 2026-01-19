using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using MathNet.Numerics.LinearAlgebra;

namespace MathTools.Core
{
    /// <summary>
    /// 主成分分析（PCA）实现，使用 SVD 进行分解。
    /// </summary>
    public sealed class BCPca
    {
        private readonly List<double[]> _patterns = new();
        private Matrix<double> _v;
        private Vector<double> _w;
        private Vector<double> _mean;
        private Vector<double> _std;

        public void AddPattern(double[] vector)
        {
            if (vector == null)
            {
                throw new ArgumentNullException(nameof(vector));
            }

            _patterns.Add((double[])vector.Clone());
        }

        public void Calculate()
        {
            if (_patterns.Count == 0)
            {
                throw new InvalidOperationException("no patterns available");
            }

            int rows = _patterns.Count;
            int cols = _patterns[0].Length;
            var matrix = Matrix<double>.Build.Dense(rows, cols);
            for (int i = 0; i < rows; i++)
            {
                if (_patterns[i].Length != cols)
                {
                    throw new InvalidOperationException("inconsistent pattern size");
                }

                for (int j = 0; j < cols; j++)
                {
                    matrix[i, j] = _patterns[i][j];
                }
            }

            _mean = matrix.ColumnSums() / rows;
            for (int i = 0; i < rows; i++)
            {
                matrix.SetRow(i, matrix.Row(i) - _mean);
            }

            _std = Vector<double>.Build.Dense(cols);
            for (int j = 0; j < cols; j++)
            {
                double variance = matrix.Column(j).PointwisePower(2.0).Sum() / (rows - 1);
                _std[j] = Math.Sqrt(variance);
                if (_std[j] > 0.0)
                {
                    matrix.SetColumn(j, matrix.Column(j) / _std[j]);
                }
            }

            var svd = matrix.Svd(true);
            _w = svd.S;
            _v = svd.VT.Transpose();
        }

        public double GetEigenvalue(int index)
        {
            if (_w == null)
            {
                throw new InvalidOperationException("PCA not calculated");
            }

            return _w[index] * _w[index];
        }

        public double[] Transform(double[] input, int components)
        {
            if (_v == null || _mean == null || _std == null)
            {
                throw new InvalidOperationException("PCA not calculated");
            }

            if (components <= 0 || components > _v.ColumnCount)
            {
                throw new ArgumentOutOfRangeException(nameof(components));
            }

            var vec = Vector<double>.Build.DenseOfArray(input);
            vec -= _mean;
            for (int i = 0; i < _std.Count; i++)
            {
                if (_std[i] > 0.0)
                {
                    vec[i] /= _std[i];
                }
            }

            var projected = _v.TransposeThisAndMultiply(vec);
            return projected.SubVector(0, components).ToArray();
        }

        public void Save(string filePath)
        {
            if (_v == null || _w == null || _mean == null || _std == null)
            {
                throw new InvalidOperationException("PCA not calculated");
            }

            using var writer = new StreamWriter(filePath);
            WriteVector(writer, _mean);
            WriteVector(writer, _std);
            WriteVector(writer, _w.PointwisePower(2.0));
            for (int i = 0; i < _v.RowCount; i++)
            {
                WriteVector(writer, _v.Row(i));
            }
        }

        public void Load(string filePath, int dimension)
        {
            using var reader = new StreamReader(filePath);
            _mean = ReadVector(reader, dimension);
            _std = ReadVector(reader, dimension);
            var eigen = ReadVector(reader, dimension);
            _w = eigen.PointwiseSqrt();

            var v = Matrix<double>.Build.Dense(dimension, dimension);
            for (int i = 0; i < dimension; i++)
            {
                var row = ReadVector(reader, dimension);
                v.SetRow(i, row);
            }

            _v = v;
        }

        private static void WriteVector(TextWriter writer, Vector<double> vector)
        {
            for (int i = 0; i < vector.Count; i++)
            {
                if (i > 0)
                {
                    writer.Write('\t');
                }

                writer.Write(vector[i].ToString("R", CultureInfo.InvariantCulture));
            }

            writer.WriteLine();
        }

        private static Vector<double> ReadVector(TextReader reader, int length)
        {
            string line = reader.ReadLine() ?? throw new InvalidOperationException("unexpected end of file");
            string[] parts = line.Split('\t', StringSplitOptions.RemoveEmptyEntries);
            if (parts.Length < length)
            {
                throw new InvalidOperationException("vector length mismatch");
            }

            var vec = Vector<double>.Build.Dense(length);
            for (int i = 0; i < length; i++)
            {
                vec[i] = double.Parse(parts[i], CultureInfo.InvariantCulture);
            }

            return vec;
        }
    }
}
