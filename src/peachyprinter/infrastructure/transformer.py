import numpy as np
import logging
logger = logging.getLogger('peachy')
from peachyprinter.domain.transformer import Transformer
import threading

class LinearAlgebraTransformer(Transformer):
    def __init__(self, upper_height, lower_points, upper_points):
        '''Given upper_height as the offset to upper points,
            <upper|lower>_points as a zippered array of (deflection,distance),
            where each <distance|deflection> are given as (x,y)
            
            This class creates a linear transform matrix that you can transform any size of inputs into.
            
            This transform works by creating a linear transform matrix
            from the matrix P (Positions from calibration) of size n=8:
            [x0 y0 z0]
            [x1 y1 z1] = P
            [  ....  ]
            [xn yn zn]

            And the matched matrix D (Deflections as a percentage of "max")
            [dx0 dy0]
            [dx1 dy1] = D
            [  ...  ]
            [dxn dyn]

            To solve the linear transforms to solve for the "Transform" matrix T:
            P.T=D
            with the sizes of: [8x3].[3x2]=[8x2]

            *NOTE: due to matrix magics, you can't just do a simple inverse
            * to solve the magic matrix equation of T=(P_left^-1).D
            *
            * So we take the left inverse, and I chose this cause it made the dimenisons work.
            *  - That's literally the only reason
            *  - The pseudo inverse in numpy seems to work too, so that's a nice bi-directional solution
            * 
            * Given a matrix A of size [5x3]
            * Left inverse solves the equation of (A-1).A = I  [3x3]
            * Right inverse solves the equation of A.(A-1) = I [5x5]
            * See the difference in sizes, that's important to making the math work.
            *** leftInverse(A) = inverse(A_t*A)*A_t

            Use the function deflections=transform(np.array([(x,y,z),(xn,yn,zn)])) to use the magic
            - Output size matches input size for an array of #rows=n

            The stupid inputs was because I didn't want to fix upstream. Feel free to fix that
            '''

        self._upper_height = upper_height
        self._lower_points = lower_points
        self._upper_points = upper_points
        self._transform = self._create_transform(lower_points,upper_points,upper_height)
        self._cache = {}

    def _get_centroids(self, points):
        '''Take set of points in [((x,y)...(x,y)), ... ] and return [(x,y),(...)] centroid list
            - This may also just be the det(points_matrix) '''
        centroids=[]
        for pointset in points:
            x = 0
            y = 0
            for (i, xy) in enumerate(pointset):
                x=x+xy[0]
                y=y+xy[1]
            total=i+1
            centroids.append((x/total,y/total))
        return centroids

    def _create_offsets(self, centroids):
        '''Takes the matrix of (upper;lower) centroids and creates the offset transform
           to satisfy the offset part of calibration (addition)'''
        print "YOU GAVE ME THIS CENTROID MATRIX:"
        print centroids

    def _create_transform(self,lower_points,upper_points,upper_height):
        '''Creates self._transform, you feed it xyz and defliections - this creates the linear transform matrix for use in self.transform()'''

        upper_distances = [distance for (deflection, distance) in upper_points.items()]
        upper_deflections = [deflection for (deflection, distance) in upper_points.items()]
        lower_distances = [distance for (deflection, distance) in lower_points.items()]
        lower_deflections = [deflection for (deflection, distance) in lower_points.items()]

        #TODO: CLean this up if the 4x4 matrix works
        #centroids = self._get_centroids((upper_points, lower_points))
        #xyz_centroids = np.concatenate((centroids, [(upper_height,),(0,)]), axis=1)
        #self._create_offsets(xyz_centroids) #Creates self._offset(deflections) used in transform

        #List of Tupples, Tupples contain each row, [colum_index][row_index]
        #[(0,)]*4 means 4 rows of 1 index each 4X1 matrix
        #[(1,1),(2,2)] is an array 2X2 with 1's on top and 2's on the bottom 
        #np.concatenate axis=1 is colum concatenate, axis=0 is concatenating a new row

        #Bring it all together into a long 8x4 array
        upper_3d_distances = np.concatenate((upper_distances, [(upper_height,)]*4), axis=1)
        lower_3d_distances = np.concatenate((lower_distances, [(0,)]*4), axis=1)

        distances_3d = np.concatenate((upper_3d_distances, lower_3d_distances), axis=0)
        distances_3d = np.concatenate((distances_3d, [(1,)]*8), axis=1)
        distances_3d_inv = self._left_inverse(distances_3d)

        #np's linalg's pseudo inverse does the same as left inverse but by solving the problem of least squares (or something)
        #distances_3d_inv = np.linalg.pinv(distances_3d)

        #Concatenate the height and constant to finish the output side of the array
        upper_deflections = np.concatenate((upper_deflections,[(upper_height,1)]*4), axis=1)
        lower_deflections = np.concatenate((lower_deflections,[(0,1)]*4), axis=1)
        deflections = np.concatenate((upper_deflections,lower_deflections), axis=0)
        transform = np.dot(distances_3d_inv, deflections)

        if (True):
            print "####### P matrix ############"
            print distances_3d
            print "####### P^-1 (left) matrix ############"
            print distances_3d_inv
            print "####### D matrix ############"
            print deflections
            print "####### T matrix ############"
            print transform
        return transform

    def _get_deflection_centroid(self,matrix):
        ''' gets the calibration centroid of set of coordinates each in a seperate row as a percent of max '''
        for coordinate in matrix:
            pass
        centroid = (0,0)
        return centroid


    def _left_inverse(self,matrix):
        ''' leftInverse(A) = inverse(A_t*A)*A_t'''

        transposed_matrix = np.transpose(matrix)
        sq_helper = np.dot(transposed_matrix,matrix)
        inverse_helper = np.linalg.inv(sq_helper)
        left_inverse = np.dot(inverse_helper,transposed_matrix)
        return left_inverse

    def transform(self,xyz_point):
        xyz_d=(xyz_point[0], xyz_point[1], xyz_point[2], 1)
        deflections = np.dot(xyz_d,self._transform)
        return deflections


class OneToOneTransformer(Transformer):
    def transform(self, xyz):
        x, y, z = xyz
        return [x, y, z]


'Takes Values from -1.0 to 1.0 on both axis and returns a scaled version between 0 and 1'


class TuningTransformer(Transformer):
    def __init__(self, scale=1.0):
        if scale > 0.0 and scale <= 1.0:
            self._scale = scale
        else:
            logger.error('Scale must be between 0.0 and 1.0 was %s' % scale)
            raise Exception('Scale must be between 0.0 and 1.0 was %s' % scale)

    def transform(self, xyz):
        x, y, z = [self._check_and_adjust(value) for value in xyz]
        x = self._transform(x)
        y = self._transform(y)
        return [x, y]

    def set_scale(self, new_scale):
        self._scale = new_scale

    def _check_and_adjust(self, value):
        if value > 1.0:
            value = 1.0
            logger.info("Adjusting Values")
        if value < 0.0:
            value = 0.0
            logger.info("Adjusting Values")
        return value

    def _transform(self, axis):
        return ((axis - 0.5) * self._scale) + 0.5


class HomogenousTransformer(Transformer):
    def __init__(self, scale, upper_height, lower_points, upper_points):
        self._lock = threading.Lock()
        self._scale = scale
        self._upper_height = upper_height

        lower_points = self._sort_points(lower_points)
        upper_points = self._sort_points(upper_points)

        deflection_scale = (
            (upper_points[0][0][0] - upper_points[2][0][0]) / (lower_points[0][0][0] - lower_points[2][0][0]),
            (upper_points[0][0][1] - upper_points[2][0][1]) / (lower_points[0][0][1] - lower_points[2][0][1])
            )

        self._lower_points = lower_points
        self._upper_points = [(self._scale_point(deflection, deflection_scale), distance) for (deflection, distance) in lower_points]

        self._get_transforms()
        self._cache = {}

    def _scale_point(self, point, scale):
        x, y = point
        dx = x * 2.0 - 1.0
        dy = y * 2.0 - 1.0
        rx = dx * scale[0]
        ry = dy * scale[1]
        tx = (rx + 1.0) / 2.0
        ty = (ry + 1.0) / 2.0
        return (tx, ty)

    def _sort_points(self, points):
        distances = [distance for (deflection, distance) in points.items()]
        normals = [(int(x / abs(x)), int(y / abs(y))) for (x, y) in distances]
        return [
            points.items()[normals.index(( 1,  1))],
            points.items()[normals.index(( 1, -1))],
            points.items()[normals.index((-1, -1))],
            points.items()[normals.index((-1,  1))],
            ]

    def _get_transforms(self):
        self._lock.acquire()
        try:
            self._lower_transform = self._get_transformation_matrix(self._lower_points)
            self._upper_transform = self._get_transformation_matrix(self._upper_points)
        finally:
            self._lock.release()

    def _get_transformation_matrix(self,mappings):
        mapping_matrix = self._build_matrix(mappings)
        b = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
        solved_matrix = np.linalg.solve(mapping_matrix,b)
        forwards = np.matrix([solved_matrix[0:3], solved_matrix[3:6], solved_matrix[6:9]])
        inverse = forwards.I
        return inverse

    def _build_matrix(self, points):
        builder = []
        index = 0
        for ((xp, yp), (xi, yi)) in points:
            augment = self._augment(index, xi / self._scale, yi / self._scale)
            builder.append([ xp, yp,  1,  0,  0,  0,  0,  0,  0] + augment[0])
            builder.append([  0,  0,  0, xp, yp,  1,  0,  0,  0] + augment[1])
            builder.append([  0,  0,  0,  0,  0,  0, xp, yp,  1] + augment[2])
            index += 1
        builder.append([  1,  1,  1,  1,  1,  1,  1,  1,  1,     0,   0,   0,   0])
        return np.array(builder)

    def _augment(self, index, xi, yi):
        augment = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
        for i in range(0, 4):
            if i == index:
                augment[0][i] = -xi
                augment[1][i] = -yi
                augment[2][i] = -1
        return augment

    def _transforms_for_height(self, height):
        if height == 0:
            return self._lower_transform
        elif height == self._upper_height:
            return self._upper_transform
        elif height in self._cache.keys():
            return self._cache[height]
        else:
            current = self._scale_transform(height)
            self._cache = {height: current}
            return current

    def _scale_transform(self, height):
        adjusted_height = height / self._upper_height
        return (adjusted_height * (self._upper_transform - self._lower_transform)) + self._lower_transform

    def transform(self, (x, y, z)):
        '''transforms cartesian x,y,z into "coneular" dimensions,
            and outputs x_power, y_power as a number between 1.0 and 0.0.
            
            This accounts for both "off center-ness" and the "cosine-ness"
            of the power curve.'''

        self._lock.acquire()
        try:
            realworld = np.array([[x], [y], [1]])
            computerland = self._transforms_for_height(z) * realworld
            [kx, ky, k] = [computerland.item(i, 0) for i in range(3)]
        finally:
            self._lock.release()
        x1, y1 = (kx/k, ky/k)
        if x1 >= 0.0 and x1 <= 1.0 and y1 >= 0.0 and y1 <= 1.0:
            return (x1, y1)
        else:
            logger.warning("Bounds of printer exceeded: %s,%s" % (x, y))
            adjusted_x = min(1.0, max(0.0, x1))
            adjusted_y = min(1.0, max(0.0, y1))
            return(adjusted_x, adjusted_y)

    def set_scale(self, new_scale):
        self._scale = new_scale
        self._get_transforms()

if __name__ == "__main__":
    #adhock debug:
    height = 50.0 

    #(deflection,distance)
    lower_points_real = { 
            (0.6807058823529412, 0.31967914438502676): (-20.0, 20.0),
            (0.6759144, 0.6752727): (20.0, 20.0),
            (0.2768556149732621, 0.6595294117647058): (20.0, -20.0),
            (0.2850695187165776, 0.3080427807486631): (-20.0, -20.0)
            }
    upper_points_real = { 
            (0.7645561497326203, 0.7457754010695187): (20.0, 20.0),
            (0.7635294117647059, 0.2402139037433155): (-20.0, 20.0),
            (0.1868449197860963, 0.22481283422459894): (-20.0, -20.0),
            (0.1680213903743315, 0.7210695187165775): (20.0, -20.0)
            }

    lower_points = { 
            (0.90, 0.90): (20.0, 20.0),
            (0.90, 0.10): (-20.0, 20.0),
            (0.10, 0.90): (20.0, -20.0),
            (0.10, 0.10): (-20.0, -20.0)
            }
    upper_points = { 
            (0.99, 0.99): (20.0, 20.0),
            (0.99, 0.01): (-20.0, 20.0),
            (0.01, 0.01): (-20.0, -20.0),
            (0.01, 0.99): (20.0, -20.0)
            }

    print "Upper points:"
    print upper_points

    print "Lower Points:"
    print lower_points
    
    example_xyz = (0.0,0.0,0.0)

    print "LinTransformerMade"
    LinTransformer=LinearAlgebraTransformer(height ,lower_points_real ,upper_points_real)
    #LinTransformer=LinearAlgebraTransformer(height ,lower_points, upper_points)
    deflections = LinTransformer.transform(example_xyz)
    print "Deflections at 0mm centered: {0}".format(deflections)

    example_xyz = (0.0,0.0,10.0)
    deflections = LinTransformer.transform(example_xyz)
    print "Deflections after 10mm centered: {0}".format(deflections)

    example_xyz = (0.0,0.0,50.0)
    deflections = LinTransformer.transform(example_xyz)
    print "Deflections after 50mm centered: {0}".format(deflections)

    import time
    now=time.time()
    iterations=100000
    for i in range(1,iterations):
        deflections = LinTransformer.transform(example_xyz)
    after=time.time()
    print "{0} iterations of transform took {1} seconds".format(iterations,after-now)

    
