using System;
using System.IO;

namespace Math.Core.Legacy.Collections
{
    /// <summary>
    /// 简化版 pair：保存两个异构值。
    /// </summary>
    public struct BCPair<TFirst, TSecond>
    {
        /// <summary>
        /// 第一个元素。
        /// </summary>
        public TFirst First;

        /// <summary>
        /// 第二个元素。
        /// </summary>
        public TSecond Second;

        /// <summary>
        /// 创建一个 pair。
        /// </summary>
        public BCPair(TFirst first, TSecond second)
        {
            First = first;
            Second = second;
        }

        /// <summary>
        /// 返回类名。
        /// </summary>
        public static string ClassName() => "BCPair<TFirst, TSecond>";

        /// <summary>
        /// 输出对象信息。
        /// </summary>
        public void Output(TextWriter writer)
        {
            if (writer == null)
            {
                throw new ArgumentNullException(nameof(writer));
            }

            writer.WriteLine($"{ClassName()} ({First}, {Second})");
        }
    }
}
