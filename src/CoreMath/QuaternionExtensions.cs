using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace CoreMath
{
    public static class QuaternionExtensions
    {
        public static float[] QuaternionFromRotationMatrix(this float[] matrix)
        {
            // http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
            // assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

            /*var m11 = matrix[0];
            var m12 = matrix[4];
            var m13 = matrix[8];
            var m21 = matrix[1];
            var m22 = matrix[5];
            var m23 = matrix[9];
            var m31 = matrix[2];
            var m32 = matrix[6];
            var m33 = matrix[10];*/
            var m11 = matrix[0]; var m12 = matrix[4]; var m13 = matrix[8];
            var m21 = matrix[1]; var m22 = matrix[5]; var m23 = matrix[9];
            var m31 = matrix[2]; var m32 = matrix[6]; var m33 = matrix[10];

            var trace = m11 + m22 + m33;
            var s = 0f;

            if (trace > 0)
            {
                s = 0.5f / (float)Math.Sqrt(trace + 1f);

                return new float[]
                {
                    (m32 - m23) * s,
                    (m13 - m31) * s,
                    (m21 - m12) * s,
                    0.25f / s
                };

            }
            else if (m11 > m22 && m11 > m33)
            {
                s = 2.0f * (float)Math.Sqrt(1f + m11 - m22 - m33);

                return new float[]
                {
                    0.25f * s,
                    (m12 + m21) / s,
                    (m13 + m31) / s,
                    (m32 - m23) / s
                };
            }
            else if (m22 > m33)
            {
                s = 2f * (float)Math.Sqrt(1.0f + m22 - m11 - m33);

                return new float[]
                {
                    (m12 + m21) / s,
                    0.25f * s,
                    (m23 + m32) / s,
                    (m13 - m31) / s
                };
            }
            else
            {
                s = 2f * (float)Math.Sqrt(1f + m33 - m11 - m22);

                return new float[]
                {
                    (m13 + m31) / s,
                    (m23 + m32) / s,
                    0.25f * s,
                    (m21 - m12) / s
                };
            }
        }


        public static float[] QuaternionProduct(this float[] a, float[] b)
        {
            return new float[]
            {
                a[0] * b[3] + a[3] * b[0] + a[1] * b[2] - a[2] * b[1],
                a[1] * b[3] + a[3] * b[1] + a[2] * b[0] - a[0] * b[2],
                a[2] * b[3] + a[3] * b[2] + a[0] * b[1] - a[1] * b[0],
                a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2]
            };
        }
        /// <summary>
        /// Creates a quaternionr from the Euler rotations specified
        /// </summary>
        /// <param name="quaternion"></param>
        /// <param name="x">X rotation in radians</param>
        /// <param name="y">Y rotation in radians</param>
        /// <param name="z">Z rotation in radians</param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static float[] QuaternionFromEuler(this float[] euler, EulerOrder order)
        {
            // http://www.mathworks.com/matlabcentral/fileexchange/
            // 	20696-function-to-convert-between-dcm-euler-angles-quaternions-and-euler-vectors/
            //	content/SpinCalc.m

            var c1 = (float)Math.Cos(euler[0] / 2);
            var c2 = (float)Math.Cos(euler[1] / 2);
            var c3 = (float)Math.Cos(euler[2] / 2);
            var s1 = (float)Math.Sin(euler[0] / 2);
            var s2 = (float)Math.Sin(euler[1] / 2);
            var s3 = (float)Math.Sin(euler[2] / 2);

            var result = new float[4];

            if (order == EulerOrder.XYZ)
            {

                result[0] = s1 * c2 * c3 + c1 * s2 * s3;
                result[1] = c1 * s2 * c3 - s1 * c2 * s3;
                result[2] = c1 * c2 * s3 + s1 * s2 * c3;
                result[3] = c1 * c2 * c3 - s1 * s2 * s3;

            }
            else if (order == EulerOrder.YXZ)
            {

                result[0] = s1 * c2 * c3 + c1 * s2 * s3;
                result[1] = c1 * s2 * c3 - s1 * c2 * s3;
                result[2] = c1 * c2 * s3 - s1 * s2 * c3;
                result[3] = c1 * c2 * c3 + s1 * s2 * s3;

            }
            else if (order == EulerOrder.ZXY)
            {

                result[0] = s1 * c2 * c3 - c1 * s2 * s3;
                result[1] = c1 * s2 * c3 + s1 * c2 * s3;
                result[2] = c1 * c2 * s3 + s1 * s2 * c3;
                result[3] = c1 * c2 * c3 - s1 * s2 * s3;

            }
            else if (order == EulerOrder.ZYX)
            {

                result[0] = s1 * c2 * c3 - c1 * s2 * s3;
                result[1] = c1 * s2 * c3 + s1 * c2 * s3;
                result[2] = c1 * c2 * s3 - s1 * s2 * c3;
                result[3] = c1 * c2 * c3 + s1 * s2 * s3;

            }
            else if (order == EulerOrder.YZX)
            {

                result[0] = s1 * c2 * c3 + c1 * s2 * s3;
                result[1] = c1 * s2 * c3 + s1 * c2 * s3;
                result[2] = c1 * c2 * s3 - s1 * s2 * c3;
                result[3] = c1 * c2 * c3 - s1 * s2 * s3;

            }
            else if (order == EulerOrder.XZY)
            {

                result[0] = s1 * c2 * c3 - c1 * s2 * s3;
                result[1] = c1 * s2 * c3 - s1 * c2 * s3;
                result[2] = c1 * c2 * s3 + s1 * s2 * c3;
                result[3] = c1 * c2 * c3 + s1 * s2 * s3;

            }

            return result;
        }


        /// <summary>
        /// Creates a new quaternion that rotates sourceDirection vector (in world space) to coincide with the
        /// targetDirection vector (in world space).
        /// Rotation is performed around the origin.
        /// The vectors sourceDirection and targetDirection are assumed to be normalized.
        /// @note There are multiple such rotations - this function returns the rotation that has the shortest angle
        /// (when decomposed to axis-angle notation).
        /// </summary>
        /// <param name="sourceDirection">Normalized original direction</param>
        /// <param name="targetDirection">Normalized target direction</param>
        /// <returns></returns>
        public static float[] QuaternionRotateFromTo(this float[] sourceDirection, float[] targetDirection)
        {
            if (sourceDirection.EqualsMatrix(targetDirection))
                return new float[] { 0, 0, 0, 1 };

            // If sourceDirection == targetDirection, the cross product comes out zero, and normalization would fail. In that case, pick an arbitrary axis.
            var axis = sourceDirection.VectorCrossProduct(targetDirection);

            float oldLength = axis.VectorLength();
            axis = axis.NormalizeVector();

            if (oldLength != 0f)
            {
                float halfCosAngle = .5f * sourceDirection.VectorDotProduct(targetDirection);
                float cosHalfAngle = (float)System.Math.Sqrt(0.5f + halfCosAngle);
                float sinHalfAngle = (float)System.Math.Sqrt(0.5f - halfCosAngle);
                return new float[] { axis[0] * sinHalfAngle, axis[1] * sinHalfAngle, axis[2] * sinHalfAngle, cosHalfAngle };
            }
            else
                return new float[] { 0f, 1f, 0f, 0f };

        }

        public static float QuaternionDotProduct(this float[] a, float[] b)
        {
            return a.X() * b.X() + a.Y() * b.Y() + a.Z() * b.Z() + a.W() * b.W();
        }

        public static float[] QuaternionSlerp(this float[] quat, float[] target, float t)
        {
            float angle = quat.QuaternionDotProduct(target);
            float sign = 1f; // Multiply by a sign of +/-1 to guarantee we rotate the shorter arc.
            if (angle < 0f)
            {
                angle = -angle;
                sign = -1f;
            }

            float a;
            float b;
            if (angle < 0.999) // perform spherical linear interpolation.
            {
                // angle = Acos(angle); // After this, angle is in the range pi/2 -> 0 as the original angle variable ranged from 0 -> 1.
                angle = (-0.69813170079773212f * angle * angle - 0.87266462599716477f) * angle + 1.5707963267948966f;

                float ta = t * angle;
                /*
                // If Sin() is based on a lookup table, prefer that over polynomial approximation.
                a = Sin(angle - ta);
                b = Sin(ta);
*/
                // Not using a lookup table, manually compute the two sines by using a very rough approximation.
                float ta2 = ta * ta;
                b = ((5.64311797634681035370e-03f * ta2 - 1.55271410633428644799e-01f) * ta2 + 9.87862135574673806965e-01f) * ta;
                a = angle - ta;
                float a2 = a * a;
                a = ((5.64311797634681035370e-03f * a2 - 1.55271410633428644799e-01f) * a2 + 9.87862135574673806965e-01f) * a;

            }
            else // If angle is close to taking the denominator to zero, resort to linear interpolation (and normalization).
            {
                a = 1f - t;
                b = t;
            }

            // Lerp and renormalize.
            return quat.VectorScale(a * sign).Add(target.VectorScale(b)).NormalizeVector();
        }

        /// <summary>
        /// Computes a quaternion rotating around a given normalized axis an angle of given radians
        /// </summary>
        /// <param name="axis">Normalized 3D axis of rotation</param>
        /// <param name="angle">Angle in radians</param>
        /// <returns></returns>
        public static float[] QuaternionFromRotation(this float[] axis, float angle)
        {
            var semiangle = angle / 2;
            var sin = (float)Math.Sin(semiangle);
            return new float[]
            {
                axis[0] * sin,
                axis[1] * sin,
                axis[2] * sin,
                (float)Math.Cos(semiangle)
            };
        }

        public static float QuaternionLength(this float[] quaternion) => (float)Math.Sqrt((quaternion[3] * quaternion[3]) + quaternion.Take(3).ToArray().VectorLengthSq());

        public static float[] NormalizeQuaternion(this float[] quaternion) => quaternion.VectorScale(1 / quaternion.QuaternionLength());

        public static float[] QuaternionToAxisAngle(this float[] q)
        {
            if (Math.Abs(q[3]) > 1.0f)
            {
               q = q.NormalizeQuaternion();
            }

            var w = 2.0f * (float)Math.Acos(q[3]);

            var den = (float)Math.Sqrt(1.0 - (q[3] * q[3]));
            if (den > 0.0001f)
            {
                return new float[]
                {
                    q[0] / den,
                    q[1] / den,
                    q[2] / den,
                    w
                };
            }
            else
            {
                // This occurs when the angle is zero.
                // Not a problem: just set an arbitrary normalized axis.
                return new float[] { 1, 0, 0, w };
            }
        }

    }
}
