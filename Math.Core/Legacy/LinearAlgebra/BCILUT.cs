using System;
using System.IO;
using MathNet.Numerics.LinearAlgebra;

namespace Math.Core.Legacy.LinearAlgebra
{
    /// <summary>
    /// ILUT 预条件器的简化实现（基于 LU 分解）。
    /// </summary>
    public sealed class BCILUT
    {
        private Matrix<double> _lower;
        private Matrix<double> _upper;
        private int _dropped;
        private double _eps;
        private int _k;

        /// <summary>
        /// 构造并初始化分解。
        /// </summary>
        public BCILUT(Matrix<double> matrix, int k, double eps)
        {
            if (matrix == null)
            {
                throw new ArgumentNullException(nameof(matrix));
            }

            _k = k;
            _eps = eps;
            Init(matrix);
        }

        /// <summary>
        /// 被丢弃元素数量（简化实现固定为 0）。
        /// </summary>
        public int Dropped => _dropped;

        /// <summary>
        /// 下三角矩阵。
        /// </summary>
        public Matrix<double> Lower => _lower;

        /// <summary>
        /// 上三角矩阵。
        /// </summary>
        public Matrix<double> Upper => _upper;

        /// <summary>
        /// 通过 LU 分解初始化。
        /// </summary>
        public void Init(Matrix<double> matrix)
        {
            var lu = matrix.LU();
            _lower = lu.L;
            _upper = lu.U;
            _dropped = 0;
        }

        /// <summary>
        /// 求解线性系统（使用 LU）。
        /// </summary>
        public Vector<double> Solve(Vector<double> rhs)
        {
            if (_lower == null || _upper == null)
            {
                throw new InvalidOperationException("ILUT not initialized");
            }

            var y = _lower.Solve(rhs);
            return _upper.Solve(y);
        }

        /// <summary>
        /// 输出对象信息。
        /// </summary>
        public void Output(TextWriter writer)
        {
            if (writer == null)
            {
                throw new ArgumentNullException(nameof(writer));
            }

            writer.WriteLine("BCILUT");
            writer.WriteLine($"Dropped: {_dropped}");
            writer.WriteLine($"Threshold: {_eps}");
            writer.WriteLine($"Fill depth: {_k}");
        }
    }
}
