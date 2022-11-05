"""
.. currentmodule:: sire.maths


"""

from ..legacy import Maths as _Maths

from .. import use_new_api as _use_new_api
_use_new_api()

Matrix = _Maths.Matrix
Quaternion = _Maths.Quaternion
Triangle = _Maths.Triangle
Torsion = _Maths.Torsion

pi = _Maths.pi

from ._vector import Vector
from ._sphere import Sphere


def create_quaternion(angle=None, axis=None,
                      matrix=None, quaternion=None):
    """Create a quaternion from the passed angle and axis
       of the passed rotation matrix. If a rotation
       matrix is passed then this will ignore the
       passed angle and axis. If a quaternion is passed
       then this will ignore the matrix, angle and axis
       arguments.

        angle: (float or angle)
            The angle to rotate by - this is interpreted as
            degrees if you pass in a float. Otherwise use
            sire.units.degrees or sire.units.radians to specify
            the angle unit. This is superseded by the
            matrix and quaternion arguments.

        axis: sire.maths.Vector (or anything that can convert to a Vector)
            The vector about which to rotate. If this is not
            specified, and no other rotation specification is
            used, then the rotation is about the z axis.
            This is superseded by the matrix and
            quaternion arguments.

        quaternion: sire.maths.Quaternion
            The Quaternion description of the rotation. Note that,
            if you pass this, then the angle, axis and matrix
            arguments will be ignored.

        matrix: sire.maths.Matrix
            The 3x3 rotation matrix that describes the rotation.
            Note that, if you pass this, then the angle and axis
            arguments will be ignored. This is superseded by
            the quaternion argument.

        Returns: sire.maths.Quaternion
            The quaternion that represents the rotation
    """
    if quaternion is None:
        if type(angle) is Quaternion:
            # the user has passed in a quaternion as the first argument
            return angle

        if type(angle is Matrix) and matrix is None:
            # the user has passed in a rotation matrix as the first argument
            matrix = angle
            angle = None
            axis = None

        if matrix is None:
            if angle is None:
                raise ValueError(
                    "You must specify either the angle, rotation matrix "
                    "or quaternion used to rotate the molecule.")

            from ..units import degrees

            try:
                angle = float(angle) * degrees
            except TypeError:
                pass

            try:
                valid_angle = angle.has_same_units(degrees)
            except Exception:
                valid_angle = False

            if not valid_angle:
                raise TypeError(
                    f"The passed angle of rotation ({angle}) has the wrong "
                    f"type ({type(angle)}). It should be an angle or a float."
                )

            if axis is None:
                axis = Vector(0, 0, 1)

            # construct from the passed angle and vector
            return Quaternion(angle, Vector(axis))
        else:
            if angle is not None or axis is not None:
                from ..utils import Console
                Console.warning(
                    "The angle and/or axis of rotation will be ignored "
                    "because you have passed in a rotation matrix.")

            from ..maths import Matrix

            if type(matrix) is not Matrix:
                raise TypeError(
                    f"The rotation matrix ({matrix}) must be of type "
                    f"sire.maths.Matrix. Type {type(matrix)} is not "
                    "supported."
                )

            return Quaternion(matrix)
    else:
        if matrix is not None:
            from ..utils import Console
            Console.warning(
                "The rotation matrix will be ignored "
                "because you have passed in a quaternion.")

        if angle is not None or axis is not None:
            from ..utils import Console
            Console.warning(
                "The angle and/or axis of rotation will be ignored "
                "because you have passed in a quaternion.")

        if type(quaternion) is not Quaternion:
            raise TypeError(
                f"The quaternion ({quaternion}) must be of type "
                f"sire.maths.Quaternion. Type {type(quaternion)} is not "
                "supported."
            )

        return quaternion
