using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Text;

namespace CoreMath.Test
{
    [TestClass]
    public class QuaternionTests
    {
        [TestMethod]
        public void TestQuaternionSlerp()
        {
            var origin = new float[] { 0, 0, 1 };
            var current = origin.QuaternionRotateFromTo(new float[] { 1, 0, 1 }.NormalizeVector());
            var final = new float[] { 1, 0, 0 };
            
            var quat = origin.QuaternionRotateFromTo(final);

            var rotated = origin.VectorTransform(quat.MatrixFromQuaternion());

            Assert.IsTrue(final.EqualsMatrix(rotated, .0001f));

            var slerp = current.QuaternionSlerp(quat, .5f);
        }
    }
}
