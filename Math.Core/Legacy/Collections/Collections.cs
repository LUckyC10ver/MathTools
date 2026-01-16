using System;
using System.IO;

namespace Math.Core.Legacy.Collections
{
    /// <summary>
    /// 简化版 pair：保存两个异构值。
    /// </summary>
    public struct BCPair<TFirst, TSecond>
    {
        public TFirst First;
        public TSecond Second;

        public BCPair(TFirst first, TSecond second)
        {
            First = first;
            Second = second;
        }

        public static string ClassName() => "BCPair<TFirst, TSecond>";

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
