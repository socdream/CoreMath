using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace CoreMath
{
    public static class VectorExtensions
    {
        public static float X(this float[] vector)
        {
            return vector[0];
        }
        public static float Y(this float[] vector)
        {
            return vector[1];
        }
        public static float Z(this float[] vector)
        {
            return vector[2];
        }
        public static float W(this float[] vector)
        {
            return vector[3];
        }

        public static int X(this int[] vector)
        {
            return vector[0];
        }
        public static int Y(this int[] vector)
        {
            return vector[1];
        }
        public static int Z(this int[] vector)
        {
            return vector[2];
        }
        public static int W(this int[] vector)
        {
            return vector[3];
        }

        /// <summary>
        /// Performs the cross product between 2 3d vectors
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static float[] VectorCrossProduct(this float[] a, float[] b)
        {
            return new float[]
            {
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0]
            };
        }

        /// <summary>
        /// Computes the dot product of this and the given vector.
        /// The dot product has a geometric interpretation of measuring how close two direction vectors are to pointing
        /// in the same direction, computing angles between vectors, or the length of a projection of one vector to another.
        /// </summary>
        /// <param name="v"></param>
        /// <returns>x*v.x + y*v.y + z*v.z</returns>
        public static float VectorDotProduct(this float[] a, float[] b)
        {
            return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        }

        /// <summary>
        /// normalizes a vector of any length
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static float[] NormalizeVector(this float[] vector)
        {
            var length = vector.VectorLength();
            var result = new float[vector.Length];

            for (int i = 0; i < vector.Length; i++)
                result[i] = vector[i] / length;

            return result;
        }

        /// <summary>
        /// Gets the squared length of a vector of any length
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static float VectorLengthSq(this float[] vector)
        {
            var res = 0f;

            for (var i = 0; i < vector.Length; i++)
                res += vector[i] * vector[i];

            return res;
        }

        /// <summary>
        /// Gets the length of a vector of any length
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static float VectorLength(this float[] vector)
        {
            return (float)Math.Sqrt(vector.VectorLengthSq());
        }

        /// <summary>
        /// Returns if a vector is normalized
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="epsilon"></param>
        /// <returns></returns>
        public static bool IsVectorNormalized(this float[] vector, float epsilon = float.Epsilon)
        {
            return vector.VectorLengthSq() < epsilon;
        }

        /// <summary>
        /// Transform the given vertex using a 4x4 matrix
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="matrix"></param>
        /// <returns></returns>
        public static float[] VectorTransform(this float[] vector, float[] matrix)
        {
            return new float[] {
                matrix[0] * vector[0] + matrix[4] * vector[1] + matrix[8]  * vector[2] + matrix[12],
                matrix[1] * vector[0] + matrix[5] * vector[1] + matrix[9] * vector[2] + matrix[13],
                matrix[2] * vector[0] + matrix[6] * vector[1] + matrix[10] * vector[2] + matrix[14]
            };
        }

        /// <summary>
        /// Extract translation from a 4x4 matrix
        /// </summary>
        /// <param name="matrix"></param>
        /// <returns></returns>
        public static float[] TranslationFromMatrix(this float[] matrix)
        {
            return new float[]
            {
                matrix[12],
                matrix[13],
                matrix[14]
            };
        }

        /// <summary>
        /// Extract the scaling from a 4x4 matrix
        /// </summary>
        /// <param name="matrix"></param>
        /// <returns></returns>
        public static float[] ScalingFromMatrix(this float[] matrix)
        {
            var sx = new float[] { matrix[0], matrix[1], matrix[2] }.VectorLength();
            var sy = new float[] { matrix[4], matrix[5], matrix[6] }.VectorLength();
            var sz = new float[] { matrix[8], matrix[9], matrix[10] }.VectorLength();

            // if determine is negative, we need to invert one scale
            var det = matrix.MatrixDeterminant();
            if (det < 0)
                sx = -sx;

            return new float[]
            {
                sx,
                sy,
                sz
            };
        }
        
        /// <summary>
        /// Convert degrees to radians
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        public static float[] ToRadians(this float[] value)
        {
            var result = new float[value.Length];

            for (var i = 0; i < value.Length; i++)
                result[i] = value[i].ToRadians();

            return result;
        }

        /// <summary>
        /// Convert gradians to degrees
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        public static float[] ToDegrees(this float[] value)
        {
            var result = new float[value.Length];

            for (var i = 0; i < value.Length; i++)
                result[i] = value[i].ToDegrees();

            return result;
        }

        /// <summary>
        /// Adds 2 vectors of the same length
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static float[] Add(this float[] a, float[] b)
        {
            var result = new float[a.Length];

            for (var i = 0; i < a.Length; i++)
                result[i] = a[i] + b[i];

            return result;
        }
        public static float[] MultiplyComponents(this float[] a, float[] b)
        {
            var result = new float[a.Length];

            for (var i = 0; i < a.Length; i++)
                result[i] = a[i] * b[i];

            return result;
        }

        /// <summary>
        /// Substracts 2 vectors of the same length
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static float[] Substract(this float[] a, float[] b)
        {
            var result = new float[a.Length];

            for (var i = 0; i < a.Length; i++)
                result[i] = a[i] - b[i];

            return result;
        }

        /// <summary>
        /// Scales a vector by the given vector (vector * factor)
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="factor"></param>
        /// <returns></returns>
        public static float[] VectorScale(this float[] vector, float factor)
        {
            var result = new float[vector.Length];

            for (var i = 0; i < vector.Length; i++)
                result[i] = vector[i] * factor;

            return result;
        }

        /// <summary>
        /// Extracts the angle in radians between 2 vectors
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static float AngleBetweenVectors(this float[] a, float[] b)
        {
            var cosa = a.VectorDotProduct(b) / (float)Math.Sqrt(a.VectorLengthSq() * b.VectorLengthSq());
            if (cosa >= 1f)
                return 0f;
            else if (cosa <= -1f)
                return (float)Math.PI;
            else
                return (float)Math.Acos(cosa);
        }

        /// <summary>
        /// Extracts the angle in radians between 2 normalized vectors
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static float AngleBetweenVectorsNormalized(this float[] a, float[] b)
        {
            var cosa = a.VectorDotProduct(b);

            if (cosa >= 1f)
                return 0f;

            else if (cosa <= -1f)
                return (float)Math.PI;
            else
                return (float)Math.Acos(cosa);
        }

        /// <summary>
        /// Returns an axis aligned with the given vector
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static List<float[]> VectorAxis(this float[] vector)
        {
            var v0 = vector.NormalizeVector();

            var v1 = new float[0];

            if ((v0.Z() < -0.01) || (v0.Z() > 0.01))
                v1 = new float[] { v0.Z(), v0.Z(), -v0.X() - v0.Y() }.NormalizeVector();
            else
                v1 = new float[] { -v0.Y() - v0.Z(), v0.X(), v0.X() }.NormalizeVector();

            var v2 = v1.VectorCrossProduct(v0).NormalizeVector();

            return new List<float[]> { v0, v1, v2 };
        }
    }
}
