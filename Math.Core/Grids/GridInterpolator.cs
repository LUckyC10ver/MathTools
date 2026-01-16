using System;
using System.Collections.Generic;
using System.IO;

namespace MathTools.Core
{
    /// <summary>
    /// 多维插值器基类，负责范围检查与输入裁剪。
    /// </summary>
    public class BCGridInterpolator
    {
        private readonly BCGridBase _grid;
        private readonly double[] _values;

        /// <summary>
        /// 使用指定网格创建插值器。
        /// </summary>
        public BCGridInterpolator(BCGridBase grid)
        {
            _grid = grid ?? throw new ArgumentNullException(nameof(grid));
            if (_grid.Dimension < 0)
            {
                throw new InvalidOperationException("grid initialized in wrong dimension");
            }

            _values = new double[_grid.Dimension];
        }

        /// <summary>
        /// 返回类名称。
        /// </summary>
        public virtual string ClassName => "BCGridInterpolator";

        /// <summary>
        /// 输出对象结构信息。
        /// </summary>
        public virtual void Output(TextWriter writer)
        {
            if (writer == null)
            {
                throw new ArgumentNullException(nameof(writer));
            }

            writer.WriteLine(ClassName);
            writer.WriteLine($"Dimension: {_grid.Dimension}");
            writer.WriteLine($"Size: {_grid.Size}");
        }

        /// <summary>
        /// 通过索引访问网格值。
        /// </summary>
        public double this[int index] => _grid.At(index);

        /// <summary>
        /// 网格点总数。
        /// </summary>
        public int Size => _grid.Size;

        /// <summary>
        /// 指定维度的最大坐标。
        /// </summary>
        public double XMax(int dimension) => _grid.XMax(dimension);

        /// <summary>
        /// 指定维度的最小坐标。
        /// </summary>
        public double XMin(int dimension) => _grid.XMin(dimension);

        /// <summary>
        /// 获取裁剪后的输入向量，供派生类使用。
        /// </summary>
        protected IReadOnlyList<double> Values => _values;

        /// <summary>
        /// 进行范围检查并裁剪输入值。
        /// </summary>
        public virtual double Evaluate(IReadOnlyList<double> input)
        {
            if (input == null)
            {
                throw new ArgumentNullException(nameof(input));
            }

            if (input.Count != _grid.Dimension)
            {
                throw new InvalidOperationException("dimension of input differs from spatial dimension");
            }

            for (int dim = 0; dim < _grid.Dimension; dim++)
            {
                double value = input[dim];
                double xmin = _grid.XMin(dim);
                double xmax = _grid.XMax(dim);
                if (value < xmin)
                {
                    _values[dim] = xmin;
                }
                else if (value >= xmax)
                {
                    _values[dim] = (1.0 - Math.Sign(xmax) * 1e-4) * xmax;
                }
                else
                {
                    _values[dim] = value;
                }
            }

            return 0.0;
        }
    }
}
