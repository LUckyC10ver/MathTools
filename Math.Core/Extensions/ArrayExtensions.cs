namespace MathTools.Core
{
    public static class ArrayExtensions
    {
        /// <summary>
        /// Initialize a one-dimensional array.
        /// </summary>
        public static T[] Resize<T>(this T[] vec, int size, T val = default)
        {
            vec = new T[size];
            for (var i = 0; i < vec.Length; i++)
            {
                vec[i] = val;
            }

            return vec;
        }

        /// <summary>
        /// Initialize a two-dimensional jagged array.
        /// </summary>
        public static T[][] Resize<T>(this T[][] vec, int nr, int nc, T val = default)
        {
            vec = new T[nr][];
            for (var i = 0; i < nr; i++)
            {
                vec[i] = new T[nc];
                for (var j = 0; j < nc; j++)
                {
                    vec[i][j] = val;
                }
            }

            return vec;
        }

        /// <summary>
        /// Initialize a three-dimensional jagged array.
        /// </summary>
        public static T[][][] Resize<T>(this T[][][] vec, int nr, int nc, int nh, T val = default)
        {
            vec = new T[nr][][];
            for (var i = 0; i < nr; i++)
            {
                vec[i] = new T[nc][];
                for (var j = 0; j < nc; j++)
                {
                    vec[i][j] = new T[nh];
                    for (var k = 0; k < nh; k++)
                    {
                        vec[i][j][k] = val;
                    }
                }
            }

            return vec;
        }

        /// <summary>
        /// Initialize a four-dimensional jagged array.
        /// </summary>
        public static T[][][][] Resize<T>(this T[][][][] vec, int nr, int nc, int nh, int nd, T val = default)
        {
            vec = new T[nr][][][];
            for (var i = 0; i < nr; i++)
            {
                vec[i] = new T[nc][][];
                for (var j = 0; j < nc; j++)
                {
                    vec[i][j] = new T[nh][];
                    for (var k = 0; k < nh; k++)
                    {
                        vec[i][j][k] = new T[nd];
                        for (var m = 0; m < nd; m++)
                        {
                            vec[i][j][k][m] = val;
                        }
                    }
                }
            }

            return vec;
        }

        /// <summary>
        /// Set all values in a one-dimensional array to a value.
        /// </summary>
        public static void Set<T>(this T[] vec, T val)
        {
            for (var i = 0; i < vec.Length; i++)
            {
                vec[i] = val;
            }
        }

        /// <summary>
        /// Set all values in a two-dimensional jagged array to a value.
        /// </summary>
        public static void Set<T>(this T[][] mat, T val)
        {
            for (var i = 0; i < mat.Length; i++)
            {
                for (var j = 0; j < mat[i].Length; j++)
                {
                    mat[i][j] = val;
                }
            }
        }
    }
}
