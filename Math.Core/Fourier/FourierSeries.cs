using System;
using System.Numerics;

namespace MathTools.Core.Fourier
{
    /// <summary>
    /// 傅里叶系数基类，提供系数存储与基本操作。
    /// </summary>
    public class FourierSeries
    {
        protected Complex[] m_data;

        /// <summary>
        /// 默认构造函数，生成空系数序列。
        /// </summary>
        public FourierSeries()
        {
            m_data = Array.Empty<Complex>();
        }

        /// <summary>
        /// 构造指定长度的系数序列。
        /// </summary>
        /// <param name="count">系数数量。</param>
        /// <param name="initValue">初始化值。</param>
        public FourierSeries(int count, Complex initValue)
        {
            m_data = new Complex[count];
            for (int i = 0; i < count; i++)
            {
                m_data[i] = initValue;
            }
        }

        /// <summary>
        /// 系数数量。
        /// </summary>
        public int Size => m_data.Length;

        /// <summary>
        /// 读取指定索引处的系数。
        /// </summary>
        public Complex this[int index]
        {
            get => m_data[index];
            set => m_data[index] = value;
        }

        /// <summary>
        /// 清空系数。
        /// </summary>
        public void Clear()
        {
            m_data = Array.Empty<Complex>();
        }

        /// <summary>
        /// 根据系数计算函数值。
        /// </summary>
        public double GetValue(double dx)
        {
            if (m_data.Length == 0)
            {
                return 0.0;
            }

            double result = m_data[0].Real;
            Complex i = new Complex(0.0, 1.0);
            for (int index = 1; index < m_data.Length; index++)
            {
                result += 2.0 * (m_data[index] * Complex.Exp(i * (index * dx))).Real;
            }

            return result;
        }

        /// <summary>
        /// 系数镜像（复共轭）。
        /// </summary>
        public void Mirror()
        {
            for (int col = 0; col < m_data.Length; col++)
            {
                m_data[col] = Complex.Conjugate(m_data[col]);
            }
        }

        /// <summary>
        /// 进行角度平移。
        /// </summary>
        public void Move(double dx)
        {
            Complex i = new Complex(0.0, 1.0);
            for (int n = 0; n < m_data.Length; n++)
            {
                m_data[n] *= Complex.Exp(-i * (n * dx));
            }
        }

        /// <summary>
        /// 重置系数长度。
        /// </summary>
        public void Resize(int n)
        {
            if (n < 0)
            {
                throw new Exception("size must be non-negative");
            }

            Array.Resize(ref m_data, n);
        }
    }
}
