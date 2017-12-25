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

        public static float[] NormalizeVector(this float[] vector)
        {
            var length = vector.VectorLength();

            return new float[]
            {
                vector[0] / length,
                vector[1] / length,
                vector[2] / length
            };
        }

        public static float VectorLengthSq(this float[] vector)
        {
            var res = 0f;

            for (var i = 0; i < vector.Length; i++)
                res += vector[i] * vector[i];

            return res;
        }

        public static float VectorLength(this float[] vector)
        {
            return (float)Math.Sqrt(vector.VectorLengthSq());
        }

        public static float[] VectorTransform(this float[] vector, float[] matrix)
        {
            return new float[] {
                matrix[0] * vector[0] + matrix[4] * vector[1] + matrix[8]  * vector[2] + matrix[12],
                matrix[1] * vector[0] + matrix[5] * vector[1] + matrix[9] * vector[2] + matrix[13],
                matrix[2] * vector[0] + matrix[6] * vector[1] + matrix[10] * vector[2] + matrix[14]
            };
        }

        public static float[] TranslationFromMatrix(this float[] matrix)
        {
            return new float[]
            {
                matrix[12],
                matrix[13],
                matrix[14]
            };
        }

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
        
        public static float[] ToRadians(this float[] value)
        {
            var result = new float[value.Length];

            for (var i = 0; i < value.Length; i++)
                result[i] = value[i].ToRadians();

            return result;
        }

        public static float[] ToDegrees(this float[] value)
        {
            var result = new float[value.Length];

            for (var i = 0; i < value.Length; i++)
                result[i] = value[i].ToDegrees();

            return result;
        }

        public static float[] Add(this float[] a, float[] b)
        {
            var result = new float[a.Length];

            for (var i = 0; i < a.Length; i++)
                result[i] = a[i] + b[i];

            return result;
        }

        public static float[] Substract(this float[] a, float[] b)
        {
            var result = new float[a.Length];

            for (var i = 0; i < a.Length; i++)
                result[i] = a[i] - b[i];

            return result;
        }

        public static float[] VectorScale(this float[] a, float b)
        {
            var result = new float[a.Length];

            for (var i = 0; i < a.Length; i++)
                result[i] = a[i] * b;

            return result;
        }

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
    }
}
