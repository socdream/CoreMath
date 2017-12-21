using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace CoreMath
{
    public static class MatrixExtensions
    {
        public static float[] IdentityMatrix(this float[] matrix)
        {
            return new float[] {
                1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1
            };
        }

        public static bool EqualsMatrix(this float[] matrix, float[] other, float epsilon = float.Epsilon)
        {
            if (matrix.Length != other.Length)
                return false;

            for (var i = 0; i < matrix.Length; i++)
                if (Math.Abs(matrix[i] - other[i]) > epsilon)
                    return false;

            return true;
        }

        public static float[] LookAtMatrix(this float[]matrix, float[] eyePos, float[] targetPos, float[] up)
        {
            var z = new float[] { eyePos[0] - targetPos[0], eyePos[1] - targetPos[1], eyePos[2] - targetPos[2] };
            z = z.NormalizeVector();

            var x = up.VectorCrossProduct(z);

            if(x.VectorLengthSq() == 0)
            {
                // eye and target are in the same vertical
                z[3] += 0.0001f;
                x = up.VectorCrossProduct(z);
            }

            x = x.NormalizeVector();

            var y = z.VectorCrossProduct(x).NormalizeVector();

            return new float[]
            {
                x[0], y[0], z[0], 0,
                x[1], y[1], z[1], 0,
                x[2], y[2], z[2], 0,
                -((x[0] * eyePos[0]) + (x[1] * eyePos[1]) + (x[2] * eyePos[2])),
                -((y[0] * eyePos[0]) + (y[1] * eyePos[1]) + (y[2] * eyePos[2])),
                -((z[0] * eyePos[0]) + (z[1] * eyePos[1]) + (z[2] * eyePos[2])),
                1
            };
        }

        /// <summary>
        /// Creates a perspective projection matrix.
        /// </summary>
        /// <param name="fovy">Angle of the field of view in the y direction (in radians)</param>
        /// <param name="aspect">Aspect ratio of the view (width / height)</param>
        /// <param name="zNear">Distance to the near clip plane</param>
        /// <param name="zFar">Distance to the far clip plane</param>
        /// <param name="result">A projection matrix that transforms camera space to raster space</param>
        /// <exception cref="System.ArgumentOutOfRangeException">
        /// Thrown under the following conditions:
        /// <list type="bullet">
        /// <item>fovy is zero, less than zero or larger than Math.PI</item>
        /// <item>aspect is negative or zero</item>
        /// <item>zNear is negative or zero</item>
        /// <item>zFar is negative or zero</item>
        /// <item>zNear is larger than zFar</item>
        /// </list>
        /// </exception>
        public static float[] PerspectiveFieldOfViewMatrix(this float[] matrix, float fovy, float aspect, float zNear, float zFar)
        {
            if (fovy <= 0 || fovy > Math.PI)
                throw new ArgumentOutOfRangeException("fovy");
            if (aspect <= 0)
                throw new ArgumentOutOfRangeException("aspect");
            if (zNear <= 0)
                throw new ArgumentOutOfRangeException("zNear");
            if (zFar <= 0)
                throw new ArgumentOutOfRangeException("zFar");

            float yMax = zNear * (float)System.Math.Tan(0.5f * fovy);
            float yMin = -yMax;
            float xMin = yMin * aspect;
            float xMax = yMax * aspect;

            return matrix.PerspectiveOffCenterMatrix(xMin, xMax, yMin, yMax, zNear, zFar);
        }

        /// <summary>
        /// Creates an perspective projection matrix.
        /// </summary>
        /// <param name="left">Left edge of the view frustum</param>
        /// <param name="right">Right edge of the view frustum</param>
        /// <param name="bottom">Bottom edge of the view frustum</param>
        /// <param name="top">Top edge of the view frustum</param>
        /// <param name="zNear">Distance to the near clip plane</param>
        /// <param name="zFar">Distance to the far clip plane</param>
        /// <param name="result">A projection matrix that transforms camera space to raster space</param>
        /// <exception cref="System.ArgumentOutOfRangeException">
        /// Thrown under the following conditions:
        /// <list type="bullet">
        /// <item>zNear is negative or zero</item>
        /// <item>zFar is negative or zero</item>
        /// <item>zNear is larger than zFar</item>
        /// </list>
        /// </exception>
        public static float[] PerspectiveOffCenterMatrix(this float[] matrix, float left, float right, float bottom, float top, float zNear, float zFar)
        {
            if (zNear <= 0)
                throw new ArgumentOutOfRangeException("zNear");
            if (zFar <= 0)
                throw new ArgumentOutOfRangeException("zFar");
            if (zNear >= zFar)
                throw new ArgumentOutOfRangeException("zNear");

            var x = (2.0f * zNear) / (right - left);
            var y = (2.0f * zNear) / (top - bottom);
            var a = (right + left) / (right - left);
            var b = (top + bottom) / (top - bottom);
            var c = -(zFar + zNear) / (zFar - zNear);
            var d = -(2.0f * zFar * zNear) / (zFar - zNear);

            return new float[]
            {
                x,0,0,0,
                0,y,0,0,
                a,b,c,-1,
                0,0,d,0
            };
        }

        public static float[] MatrixInverse(this float[] matrix)
        {
            // assumes determinant is not 0
            // that is, the matrix does have an inverse
            var n = 4; // cols
            var result = new float[matrix.Length]; // make a copy of matrix

            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    result[i * n + j] = matrix[i * n + j];

            // combined lower & upper
            var toggle = matrix.MatrixDecompose(out float[] lum, out int[] perm);

            var b = new float[n];

            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                    if (i == perm[j])
                        b[j] = 1f;
                    else
                        b[j] = 0f;

                var x = Helper(lum, b);

                for (int j = 0; j < n; ++j)
                    result[j * n + i] = x[j];
            }

            return result;
        }

        public static float[] Helper(float[] luMatrix, float[] b) // helper
        {
            int n = 4; // cols
            var x = new float[n];
            b.CopyTo(x, 0);

            for (int i = 1; i < n; ++i)
            {
                var sum = x[i];
                for (int j = 0; j < i; ++j)
                    sum -= luMatrix[i * n + j] * x[j];
                x[i] = sum;
            }

            x[n - 1] /= luMatrix[(n - 1) * n + (n - 1)];
            for (int i = n - 2; i >= 0; --i)
            {
                var sum = x[i];

                for (int j = i + 1; j < n; ++j)
                    sum -= luMatrix[i * n + j] * x[j];

                x[i] = sum / luMatrix[i * n + i];
            }

            return x;
        }

        /// <summary>
        /// Crout's LU decomposition for matrix determinant and inverse
        /// stores combined lower & upper in lum[][]
        /// stores row permuations into perm[]
        /// returns +1 or -1 according to even or odd number of row permutations
        /// lower gets dummy 1.0s on diagonal (0.0s above)
        /// upper gets lum values on diagonal (0.0s below)
        /// </summary>
        /// <param name="matrix">This must be a 4x4 row major matrix</param>
        /// <param name="lum"></param>
        /// <param name="perm"></param>
        /// <returns></returns>
        public static int MatrixDecompose(this float[] matrix, out float[] lum, out int[] perm)
        {
            var toggle = 1; // even (+1) or odd (-1) row permutatuions

            // columns
            int n = 4;

            // make a copy of m[][] into result lu[][]
            lum = new float[matrix.Length];

            Array.Copy(matrix, lum, matrix.Length);

            // make perm[]
            perm = new int[n];

            for (int i = 0; i < n; ++i)
                perm[i] = i;

            for (int j = 0; j < n - 1; ++j) // process by column. note n-1 
            {
                var max = Math.Abs(lum[j * n + j]);
                var piv = j;

                for (int i = j + 1; i < n; ++i) // find pivot index
                {
                    var xij = Math.Abs(lum[i * n + j]);

                    if (xij > max)
                    {
                        max = xij;
                        piv = i;
                    }
                } // i

                if (piv != j)
                {
                    // swap rows j, piv
                    for (var i = 0; i < n; i++)
                    {
                        var temp = lum[piv * n + i];
                        lum[piv * n + i] = lum[j * n + i];
                        lum[j * n + i] = temp;
                    }

                    var t = perm[piv]; // swap perm elements
                    perm[piv] = perm[j];
                    perm[j] = t;

                    toggle = -toggle;
                }

                var xjj = lum[j * n + j];

                if (xjj != 0f)
                {
                    for (int i = j + 1; i < n; ++i)
                    {
                        var xij = lum[i * n + j] / xjj;
                        lum[i * n + j] = xij;

                        for (int k = j + 1; k < n; ++k)
                            lum[i * n + k] -= xij * lum[j * n + k];
                    }
                }

            } // j

            return toggle;
        }

        public static float MatrixDeterminant(this float[] matrix)
        {
            var cols = 4;

            var toggle = MatrixDecompose(matrix, out float[] lum, out int[] perm);
            var result = (float)toggle;

            for (int i = 0; i < cols; ++i)
                result *= lum[i * cols + i];

            return result;
        }

        public static float[] MatrixProduct(this float[] matrixA, float[] matrixB)
        {
            var aRows = 4; // matrixA.Length;
            var aCols = 4; // matrixA[0].Length;
            var bRows = 4; // matrixB.Length;
            var bCols = 4; // matrixB[0].Length;

            if (aCols != bRows)
                throw new Exception("Non-conformable matrices");

            var result = new float[aRows * bCols];

            for (int i = 0; i < aRows; ++i) // each row of A
                for (int j = 0; j < bCols; ++j) // each col of B
                    for (int k = 0; k < aCols; ++k) // could use k < bRows
                        result[i * aRows + j] += matrixA[i * aRows + k] * matrixB[k * aCols + j];

            return result;
        }

        public static float[] NormalizeMatrix(this float[] matrix)
        {
            var result = (float[])matrix.Clone();

            var row1 = new float[] { matrix[0], matrix[1], matrix[2] }.NormalizeVector();
            var row2 = new float[] { matrix[4], matrix[5], matrix[6] }.NormalizeVector();
            var row3 = new float[] { matrix[8], matrix[9], matrix[10] }.NormalizeVector();

            result[0] = row1[0];
            result[1] = row1[1];
            result[2] = row1[2];
            result[4] = row2[0];
            result[5] = row2[1];
            result[6] = row2[2];
            result[8] = row3[0];
            result[9] = row3[1];
            result[10] = row3[2];

            return result;
        }

        public static float[] ExtractLower(this float[] lum)
        {
            // lower part of an LU Doolittle decomposition (dummy 1.0s on diagonal, 0.0s above)
            var n = 4; //cols
            var result = new float[lum.Length];

            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (i == j)
                        result[i * n + j] = 1f;
                    else if (i > j)
                        result[i * n + j] = lum[i * n + j];
                }
            }

            return result;
        }

        public static float[] ExtractUpper(this float[] lum)
        {
            // upper part of an LU (lu values on diagional and above, 0.0s below)
            int n = 4; //cols
            var result = new float[lum.Length];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (i <= j)
                        result[i * n + j] = lum[i * n + j];
                }
            }
            return result;
        }

        public static string ToMatrixString(this float[] matrix)
        {
            var rows = 4;
            var cols = 4;

            var result = "";

            for (var i = 0; i < rows; i++)
            {

                for (var j = 0; j < cols; j++)
                    result += matrix[i * rows + j].ToString() + ", ";

                result += "\r\n";
            }

            return result;
        }

        public static float[] TransposeMatrix(this float[] matrix)
        {
            var columns = (matrix.Length == 9) ? 3 : (matrix.Length == 16) ? 4 : 0;

            if (columns == 0)
                throw new ArgumentException("Matrix has to be a 3x3 or 4x4 matrix.");

            var result = new float[matrix.Length];

            Array.Copy(matrix, result, matrix.Length);

            for (var i = 0; i < columns; i++)
                for (var j = 0; j < columns; j++)
                    result[i * columns + j] = matrix[j * columns + i];

            return result;
        }

        public static float[] TranslationMatrix(this float[] translation)
        {
            var result = new float[] { }.IdentityMatrix();

            /*result[0 * 4 + 3] = translation[0];
            result[1 * 4 + 3] = translation[1];
            result[2 * 4 + 3] = translation[2];*/
            result[12] = translation[0];
            result[13] = translation[1];
            result[14] = translation[2];

            return result;
        }

        /// <summary>
        /// Creates a scaling matrix from a 3 component vector (x, y, z)
        /// </summary>
        /// <param name="scale"></param>
        /// <returns></returns>
        public static float[] ScalingMatrix(this float[] scale)
        {
            var result = new float[] { }.IdentityMatrix();

            result[0 * 4 + 0] = scale[0];
            result[1 * 4 + 1] = scale[1];
            result[2 * 4 + 2] = scale[2];

            return result;
        }

        public static void MatrixDecompose(this float[] matrix, out float[] position, out float[] quaternion, out float[] scale)
        {
            var sx = new float[] { matrix[0], matrix[1], matrix[2] }.VectorLength();
            var sy = new float[] { matrix[4], matrix[5], matrix[6] }.VectorLength();
            var sz = new float[] { matrix[8], matrix[9], matrix[10] }.VectorLength();

            // if determine is negative, we need to invert one scale
            var det = matrix.MatrixDeterminant();
            if (det < 0)
                sx = -sx;

            position = new float[]
            {
                matrix[12],
                matrix[13],
                matrix[14]
            };

            // scale the rotation part
            var temp = new float[16];
            Array.Copy(matrix, temp, 16);

            var invSX = 1 / sx;
            var invSY = 1 / sy;
            var invSZ = 1 / sz;

            temp[0] *= invSX;
            temp[1] *= invSX;
            temp[2] *= invSX;

            temp[4] *= invSY;
            temp[5] *= invSY;
            temp[6] *= invSY;

            temp[8] *= invSZ;
            temp[9] *= invSZ;
            temp[10] *= invSZ;

            quaternion = matrix.QuaternionFromRotationMatrix();

            scale = new float[]
            {
                sx,
                sy,
                sz
            };
        }

        public static float[] MatrixCompose(this float[] matrix, float[] position, float[] quaternion, float[] scale)
        {
            var result = quaternion.MatrixFromQuaternion();

            result = result.ScaleMatrix(scale);

            return result.SetPositionMatrix(position);
        }

        public static float[] SetPositionMatrix(this float[] matrix, float[] position)
        {
            var result = (float[])matrix.Clone();

            /*result[0 * 4 + 3] = position[0];
            result[1 * 4 + 3] = position[1];
            result[2 * 4 + 3] = position[2];*/
            result[12] = position[0];
            result[13] = position[1];
            result[14] = position[2];

            return result;
        }

        /// <summary>
        /// Scales the given matrix
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="scale"></param>
        /// <returns></returns>
        public static float[] ScaleMatrix(this float[] matrix, float[] scale)
        {
            var result = (float[])matrix.Clone();

            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 4; j++)
                    result[i * 4 + j] *= scale[i];

            return result;
        }

        public static float[] MatrixFromQuaternion(this float[] quaternion)
        {
            var x = quaternion[0];
            var y = quaternion[1];
            var z = quaternion[2];
            var w = quaternion[3];
            var x2 = x + x;
            var y2 = y + y;
            var z2 = z + z;
            var xx = x * x2;
            var xy = x * y2;
            var xz = x * z2;
            var yy = y * y2;
            var yz = y * z2;
            var zz = z * z2;
            var wx = w * x2;
            var wy = w * y2;
            var wz = w * z2;

            var result = new float[16];

            result[0] = 1 - (yy + zz);
            result[4] = xy - wz;
            result[8] = xz + wy;

            result[1] = xy + wz;
            result[5] = 1 - (xx + zz);
            result[9] = yz - wx;

            result[2] = xz - wy;
            result[6] = yz + wx;
            result[10] = 1 - (xx + yy);

            // last column
            result[3] = 0;
            result[7] = 0;
            result[11] = 0;

            // bottom row
            result[12] = 0;
            result[13] = 0;
            result[14] = 0;
            result[15] = 1;

            return result;
        }

        public static float[] MatrixFromEuler(this float[] euler, EulerOrder order)
        {
            var result = new float[16];

            var x = euler[0];
            var y = euler[1];
            var z = euler[2];
            var a = (float)Math.Cos(x);
            var b = (float)Math.Sin(x);
            var c = (float)Math.Cos(y);
            var d = (float)Math.Sin(y);
            var e = (float)Math.Cos(z);
            var f = (float)Math.Sin(z);

            if (order == EulerOrder.XYZ)
            {

                var ae = a * e;
                var af = a * f;
                var be = b * e;
                var bf = b * f;

                result[0] = c * e;
                result[4] = -c * f;
                result[8] = d;

                result[1] = af + be * d;
                result[5] = ae - bf * d;
                result[9] = -b * c;

                result[2] = bf - ae * d;
                result[6] = be + af * d;
                result[10] = a * c;

            }
            else if (order == EulerOrder.YXZ)
            {

                var ce = c * e;
                var cf = c * f;
                var de = d * e;
                var df = d * f;

                result[0] = ce + df * b;
                result[4] = de * b - cf;
                result[8] = a * d;

                result[1] = a * f;
                result[5] = a * e;
                result[9] = -b;

                result[2] = cf * b - de;
                result[6] = df + ce * b;
                result[10] = a * c;

            }
            else if (order == EulerOrder.ZXY)
            {

                var ce = c * e;
                var cf = c * f;
                var de = d * e;
                var df = d * f;

                result[0] = ce - df * b;
                result[4] = -a * f;
                result[8] = de + cf * b;

                result[1] = cf + de * b;
                result[5] = a * e;
                result[9] = df - ce * b;

                result[2] = -a * d;
                result[6] = b;
                result[10] = a * c;

            }
            else if (order == EulerOrder.ZYX)
            {

                var ae = a * e;
                var af = a * f;
                var be = b * e;
                var bf = b * f;

                result[0] = c * e;
                result[4] = be * d - af;
                result[8] = ae * d + bf;

                result[1] = c * f;
                result[5] = bf * d + ae;
                result[9] = af * d - be;

                result[2] = -d;
                result[6] = b * c;
                result[10] = a * c;

            }
            else if (order == EulerOrder.YZX)
            {

                var ac = a * c;
                var ad = a * d;
                var bc = b * c;
                var bd = b * d;

                result[0] = c * e;
                result[4] = bd - ac * f;
                result[8] = bc * f + ad;

                result[1] = f;
                result[5] = a * e;
                result[9] = -b * e;

                result[2] = -d * e;
                result[6] = ad * f + bc;
                result[10] = ac - bd * f;

            }
            else if (order == EulerOrder.XZY)
            {
                var ac = a * c;
                var ad = a * d;
                var bc = b * c;
                var bd = b * d;

                result[0] = c * e;
                result[4] = -f;
                result[8] = d * e;

                result[1] = ac * f + bd;
                result[5] = a * e;
                result[9] = ad * f - bc;

                result[2] = bc * f - ad;
                result[6] = b * e;
                result[10] = bd * f + ac;

            }

            // last column
            result[3] = 0;
            result[7] = 0;
            result[11] = 0;

            // bottom row
            result[12] = 0;
            result[13] = 0;
            result[14] = 0;
            result[15] = 1;

            return result;
        }
    }
}
