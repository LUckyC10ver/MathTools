using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;

namespace MathTools.Core.Grids
{
    /// <summary>
    /// 多维矩形网格，使用行优先顺序存储值。
    /// </summary>
    public class BCMultivariateGrid : BCGrid
    {
        private readonly List<List<double>> _coordinates = new();
        private readonly List<int> _strides = new();
        private readonly List<double> _values = new();

        /// <summary>
        /// 默认构造：创建空网格。
        /// </summary>
        public BCMultivariateGrid()
        {
        }

        /// <summary>
        /// 根据各维度点数创建网格，并用 <paramref name="initializer"/> 初始化值。
        /// </summary>
        public BCMultivariateGrid(IReadOnlyList<int> noPoints, double initializer = -MathCoreInfo.BcInfinity)
        {
            if (noPoints == null)
            {
                throw new ArgumentNullException(nameof(noPoints));
            }

            if (noPoints.Count == 0)
            {
                return;
            }

            _coordinates.Clear();
            _strides.Clear();
            _values.Clear();

            int product = 1;
            for (int i = 0; i < noPoints.Count; i++)
            {
                int n = noPoints[i];
                if (n < 2)
                {
                    throw new InvalidOperationException("invalid number of points");
                }

                product *= n;
                var coord = new List<double>(n);
                double step = 1.0 / (n - 1);
                double x = 0.0;
                for (int j = 0; j < n; j++)
                {
                    coord.Add(x);
                    x += step;
                }

                _coordinates.Add(coord);
                _strides.Add(0);
            }

            ComputeStrides();
            for (int i = 0; i < product; i++)
            {
                _values.Add(initializer);
            }
        }

        /// <inheritdoc />
        public override int Dimension => _coordinates.Count;

        /// <inheritdoc />
        public override IReadOnlyList<double> Values => _values;

        /// <inheritdoc />
        public override bool Empty => _values.Count == 0;

        /// <inheritdoc />
        public override int Size => _values.Count;

        /// <summary>
        /// 读取网格数据（包含坐标与值）。
        /// </summary>
        public void Load(TextReader reader)
        {
            var tokenReader = new TokenReader(reader);
            int size = tokenReader.ReadInt();
            int dimension = tokenReader.ReadInt();

            if (size < 0 || dimension < 0)
            {
                throw new InvalidOperationException("invalid dimensions");
            }

            _values.Clear();
            _coordinates.Clear();
            _strides.Clear();

            for (int i = 0; i < dimension; i++)
            {
                List<double> coords = tokenReader.ReadDoubleList();
                coords.Sort();
                _coordinates.Add(coords);
                _strides.Add(0);
            }

            ComputeStrides();
            List<double> values = tokenReader.ReadDoubleList();
            if (values.Count != size)
            {
                throw new InvalidOperationException("value size mismatch");
            }

            _values.AddRange(values);
        }

        /// <summary>
        /// 写入网格数据（包含坐标与值）。
        /// </summary>
        public void Store(TextWriter writer)
        {
            if (writer == null)
            {
                throw new ArgumentNullException(nameof(writer));
            }

            writer.WriteLine(Size);
            writer.WriteLine(Dimension);

            foreach (List<double> coord in _coordinates)
            {
                writer.Write(coord.Count.ToString(CultureInfo.InvariantCulture));
                for (int i = 0; i < coord.Count; i++)
                {
                    writer.Write(' ');
                    writer.Write(coord[i].ToString("R", CultureInfo.InvariantCulture));
                }

                writer.WriteLine();
            }

            writer.WriteLine(_values.Count.ToString(CultureInfo.InvariantCulture));
            foreach (double value in _values)
            {
                writer.WriteLine(value.ToString("R", CultureInfo.InvariantCulture));
            }
        }

        /// <inheritdoc />
        public override int SizeOfDimension(int dimension)
        {
            ValidateDimension(dimension);
            return _coordinates[dimension].Count;
        }

        /// <inheritdoc />
        public override double XMax(int dimension)
        {
            ValidateDimension(dimension);
            var coords = _coordinates[dimension];
            return coords[^1];
        }

        /// <inheritdoc />
        public override double XMin(int dimension)
        {
            ValidateDimension(dimension);
            var coords = _coordinates[dimension];
            return coords[0];
        }

        /// <inheritdoc />
        public override double GetCoordinate(int dimension, int index)
        {
            ValidateDimension(dimension);
            var coords = _coordinates[dimension];
            if (index < 0 || index >= coords.Count)
            {
                throw new ArgumentOutOfRangeException(nameof(index));
            }

            return coords[index];
        }

        /// <inheritdoc />
        public override IReadOnlyList<double> GetPoints(int dimension)
        {
            ValidateDimension(dimension);
            return _coordinates[dimension];
        }

        /// <inheritdoc />
        public override int GetIndex(IReadOnlyList<int> multiIndex)
        {
            if (multiIndex == null)
            {
                throw new ArgumentNullException(nameof(multiIndex));
            }

            if (multiIndex.Count != Dimension)
            {
                throw new InvalidOperationException("incompatible dimension");
            }

            int index = 0;
            for (int i = 0; i < multiIndex.Count; i++)
            {
                int component = multiIndex[i];
                if (component < 0 || component >= _coordinates[i].Count)
                {
                    throw new InvalidOperationException("index out of bounds");
                }

                index += component * _strides[i];
            }

            return index;
        }

        /// <inheritdoc />
        public override IReadOnlyList<int> GetMultiindex(int index)
        {
            if (index < 0 || index >= _values.Count)
            {
                throw new ArgumentOutOfRangeException(nameof(index));
            }

            var multiIndex = new int[Dimension];
            int remaining = index;
            for (int i = 0; i < Dimension; i++)
            {
                int stride = _strides[i];
                multiIndex[i] = remaining / stride;
                remaining %= stride;
            }

            return multiIndex;
        }

        /// <inheritdoc />
        public override IReadOnlyList<int> GetMaximalMultiindex(IReadOnlyList<double> coordinates)
        {
            if (coordinates == null)
            {
                throw new ArgumentNullException(nameof(coordinates));
            }

            if (coordinates.Count != Dimension)
            {
                throw new InvalidOperationException("wrong dimension");
            }

            var index = new int[Dimension];
            for (int i = 0; i < Dimension; i++)
            {
                index[i] = FindIntervalIndex(_coordinates[i], coordinates[i]);
            }

            return index;
        }

        /// <inheritdoc />
        public override IReadOnlyList<int> GetMultiindex(IReadOnlyList<double> coordinates, double maxError = 2e-2)
        {
            if (coordinates == null)
            {
                throw new ArgumentNullException(nameof(coordinates));
            }

            if (coordinates.Count < Dimension)
            {
                throw new InvalidOperationException("wrong dimension");
            }

            var index = new List<int>(Dimension);
            for (int i = 0; i < Dimension; i++)
            {
                int matchIndex = -1;
                double target = coordinates[i];
                var coordList = _coordinates[i];
                for (int j = 0; j < coordList.Count; j++)
                {
                    double coord = coordList[j];
                    if (coord != 0.0)
                    {
                        if (Math.Abs((coord - target) / coord) < maxError)
                        {
                            matchIndex = j;
                            break;
                        }
                    }
                    else if (Math.Abs(target) < maxError)
                    {
                        matchIndex = j;
                        break;
                    }
                }

                if (matchIndex < 0)
                {
                    return Array.Empty<int>();
                }

                index.Add(matchIndex);
            }

            return index;
        }

        /// <inheritdoc />
        public override double GetPoint(IReadOnlyList<int> multiIndex)
        {
            if (multiIndex == null)
            {
                throw new ArgumentNullException(nameof(multiIndex));
            }

            int index = 0;
            for (int i = 0; i < multiIndex.Count; i++)
            {
                index += multiIndex[i] * _strides[i];
            }

            return _values[index];
        }

        /// <inheritdoc />
        public override double MaxValue()
        {
            if (_values.Count == 0)
            {
                return -MathCoreInfo.BcInfinity;
            }

            double max = -MathCoreInfo.BcInfinity;
            foreach (double value in _values)
            {
                if (value > max)
                {
                    max = value;
                }
            }

            return max;
        }

        /// <inheritdoc />
        public override double MinValue()
        {
            if (_values.Count == 0)
            {
                return MathCoreInfo.BcInfinity;
            }

            double min = MathCoreInfo.BcInfinity;
            foreach (double value in _values)
            {
                if (value > -MathCoreInfo.BcInfinity && value < min)
                {
                    min = value;
                }
            }

            return min;
        }

        /// <inheritdoc />
        public override double MeanValue()
        {
            if (_values.Count == 0)
            {
                return -MathCoreInfo.BcInfinity;
            }

            double sum = 0.0;
            int count = 0;
            foreach (double value in _values)
            {
                if (value > -MathCoreInfo.BcInfinity && value < MathCoreInfo.BcInfinity)
                {
                    sum += value;
                    count++;
                }
            }

            if (count == 0)
            {
                return 0.0;
            }

            return sum / count;
        }

        /// <inheritdoc />
        public override void Output(TextWriter writer)
        {
            if (writer == null)
            {
                throw new ArgumentNullException(nameof(writer));
            }

            writer.WriteLine("BCMultivariateGrid");
            writer.WriteLine($"Dimension: {Dimension}");
            writer.WriteLine($"Size: {Size}");
        }

        /// <inheritdoc />
        public override void ReadCoordinates(TextReader reader)
        {
            var tokenReader = new TokenReader(reader);
            int size = tokenReader.ReadInt();
            int dimension = tokenReader.ReadInt();
            if (size < 0 || dimension < 0)
            {
                throw new InvalidOperationException("invalid dimensions");
            }

            _coordinates.Clear();
            _strides.Clear();
            for (int i = 0; i < dimension; i++)
            {
                var coords = tokenReader.ReadDoubleList();
                coords.Sort();
                _coordinates.Add(coords);
                _strides.Add(0);
            }

            int product = 1;
            foreach (List<double> coord in _coordinates)
            {
                product *= coord.Count;
            }

            if (product != size)
            {
                throw new InvalidOperationException("product of division sizes does not match value size");
            }

            ComputeStrides();
        }

        /// <inheritdoc />
        public override int ReadPoints(TextReader reader)
        {
            if (reader == null)
            {
                throw new ArgumentNullException(nameof(reader));
            }

            var tokenReader = new TokenReader(reader);
            int number = 0;
            while (true)
            {
                if (!tokenReader.TryReadDouble(out double first))
                {
                    break;
                }

                var point = new double[Dimension + 1];
                point[0] = first;
                for (int i = 1; i < point.Length; i++)
                {
                    point[i] = tokenReader.ReadDouble();
                }

                var multiIndex = GetMultiindex(point);
                if (multiIndex.Count > 0)
                {
                    int index = GetIndex(multiIndex);
                    EnsureValueSize();
                    _values[index] = point[^1];
                }

                number++;
            }

            return number;
        }

        /// <inheritdoc />
        public override void SetPoints(int dimension, IReadOnlyList<double> coord)
        {
            ValidateDimension(dimension);
            _coordinates[dimension] = new List<double>(coord);
            EnsureValueSize();
        }

        /// <summary>
        /// 调整维度数量。
        /// </summary>
        public void SetDimension(int dimension)
        {
            if (dimension < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(dimension));
            }

            _coordinates.Clear();
            _strides.Clear();
            for (int i = 0; i < dimension; i++)
            {
                _coordinates.Add(new List<double>());
                _strides.Add(0);
            }
        }

        /// <summary>
        /// 重新计算行优先步长。
        /// </summary>
        public void ComputeStrides()
        {
            if (_coordinates.Count == 0)
            {
                return;
            }

            int product = 1;
            for (int i = _coordinates.Count - 1; i >= 0; i--)
            {
                _strides[i] = product;
                product *= _coordinates[i].Count;
            }
        }

        /// <inheritdoc />
        public override double this[int index] => _values[index];

        private void EnsureValueSize()
        {
            int product = 1;
            foreach (List<double> coord in _coordinates)
            {
                product *= coord.Count;
            }

            if (_values.Count != product)
            {
                _values.Clear();
                for (int i = 0; i < product; i++)
                {
                    _values.Add(-MathCoreInfo.BcInfinity);
                }
            }
        }

        private void ValidateDimension(int dimension)
        {
            if (dimension < 0 || dimension >= Dimension)
            {
                throw new ArgumentOutOfRangeException(nameof(dimension));
            }
        }

        private static int FindIntervalIndex(IReadOnlyList<double> vector, double x)
        {
            if (vector.Count == 0)
            {
                return -1;
            }

            if (x < vector[0])
            {
                return -1;
            }

            int n = vector.Count;
            if (x > vector[n - 1])
            {
                return n;
            }

            int klo = 0;
            int khi = n;
            while (khi - klo > 1)
            {
                int k = (khi + klo) / 2;
                if (vector[k] > x)
                {
                    khi = k;
                }
                else
                {
                    klo = k;
                }
            }

            return klo;
        }

        private sealed class TokenReader
        {
            private readonly TextReader _reader;
            private readonly Queue<string> _buffer = new();

            public TokenReader(TextReader reader)
            {
                _reader = reader ?? throw new ArgumentNullException(nameof(reader));
            }

            public bool TryReadDouble(out double value)
            {
                string token = NextToken();
                if (token == null)
                {
                    value = 0.0;
                    return false;
                }

                value = double.Parse(token, CultureInfo.InvariantCulture);
                return true;
            }

            public int ReadInt()
            {
                string token = NextToken() ?? throw new InvalidOperationException("unexpected end of stream");
                return int.Parse(token, CultureInfo.InvariantCulture);
            }

            public double ReadDouble()
            {
                string token = NextToken() ?? throw new InvalidOperationException("unexpected end of stream");
                return double.Parse(token, CultureInfo.InvariantCulture);
            }

            public List<double> ReadDoubleList()
            {
                int count = ReadInt();
                if (count < 0)
                {
                    throw new InvalidOperationException("invalid vector size");
                }

                var list = new List<double>(count);
                for (int i = 0; i < count; i++)
                {
                    list.Add(ReadDouble());
                }

                return list;
            }

            private string NextToken()
            {
                while (_buffer.Count == 0)
                {
                    string line = _reader.ReadLine();
                    if (line == null)
                    {
                        return null;
                    }

                    foreach (string token in line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries))
                    {
                        _buffer.Enqueue(token);
                    }
                }

                return _buffer.Dequeue();
            }
        }
    }
}
