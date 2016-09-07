from numpy.linalg import inv
from numpy.linalg import pinv
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
            from the matrix C (cartesian coordinates from calibration) of size n=8:
            [x0 y0 z0]
            [x1 y1 z1] = C
            [  ....  ]
            [xn yn zn]

            And the matched matrix D (Deflections as a percentage of "max")
            [dx0 dy0]
            [dx1 dy1] = D
            [  ...  ]
            [dxn dyn]

            Z matrix used in offseting laser position
            [z0 1]
            [z1 1] = Z
            [ .. ]
            [zn 1]            

            To solve the linear transforms to solve for the "Transform" matrix G
            D = C.G + Z.M
            -where the M matrix is derived from the centroids of top and bottom cartesian.

            *NOTE: due to matrix magics, you can't just do a simple inverse
            * to solve the magic matrix equation of T=(P_left^-1).D
            *
            * So we take the left inverse, this gives us the correct dimensions
            *  - the moore-penrose approach seems to be working?
            *  - The pseudo inverse in numpy seems to work too.
            * 
            *** leftInverse(A) = inverse(A_t*A)*A_t

            Use the function deflections=transform(np.array([(x,y,z),(xn,yn,zn)])) to use the magic
            - Output size matches input size for an array of #rows=n

            The stupid inputs was because I didn't want to fix upstream. Feel free to fix that
            '''

        self._upper_height = upper_height
        self._lower_points = lower_points
        self._upper_points = upper_points
        self._transform = self._create_transform(lower_points,upper_points,upper_height)

    def _create_transform(self,lower_points,upper_points,upper_height):
        '''Creates self._transform, self._cartesian_offsets, and self._deflection_offsets '''

        deflection_averages=[0,0,0]
        deflection_matrix = []
        deflection_normalized_matrix = []
        cartesian_averages=[0,0,0]
        cartesian_matrix = []
        cartesian_normalized_matrix = []
        lower_height = 0.0
        num_points = 8.0

        for (deflections, distances) in upper_points.items():
            deflection_averages[0] += deflections[0]/num_points
            deflection_averages[1] += deflections[1]/num_points
            deflection_averages[2] += upper_height/num_points
            deflection_matrix.append([deflections[0], deflections[1], upper_height])
            cartesian_averages[0] += distances[0]/num_points
            cartesian_averages[1] += distances[1]/num_points
            cartesian_averages[2] += upper_height/num_points
            cartesian_matrix.append([distances[0], distances[1], upper_height])

        for (deflections, distances) in lower_points.items():
            deflection_averages[0] += deflections[0]/num_points
            deflection_averages[1] += deflections[1]/num_points
            deflection_averages[2] += lower_height/num_points
            deflection_matrix.append([deflections[0], deflections[1], lower_height])
            cartesian_averages[0] += distances[0]/num_points
            cartesian_averages[1] += distances[1]/num_points
            cartesian_averages[2] += lower_height/num_points
            cartesian_matrix.append([distances[0], distances[1], lower_height])

        deflection_offset_matrix = np.matrix([(deflection_averages[0], deflection_averages[1], deflection_averages[2])]*8)
        cartesian_offset_matrix = np.matrix([(cartesian_averages[0], cartesian_averages[1], cartesian_averages[2])]*8)

        deflection_normalized_matrix = deflection_matrix - deflection_offset_matrix
        cartesian_normalized_matrix = cartesian_matrix - cartesian_offset_matrix

        deflections_squared = np.square(deflection_normalized_matrix)
        cartesians_squared = np.square(cartesian_normalized_matrix)

        cartesian_vandermonde = np.concatenate((cartesian_normalized_matrix, cartesians_squared), axis=1)
        cartesian_vandermonde = np.concatenate((cartesian_vandermonde, [(1,)]*8), axis=1)
        deflection_vandermonde = np.concatenate((deflection_normalized_matrix, deflections_squared), axis=1)
        deflection_vandermonde = np.concatenate((deflection_vandermonde, [(1,)]*8), axis=1)

        print cartesian_vandermonde
        cartesian_vandermonde = np.delete(cartesian_vandermonde,7,0)
        #cartesian_vandermonde = np.delete(cartesian_vandermonde,3,0)
        deflection_vandermonde = np.delete(deflection_vandermonde,7,0)
        #deflection_vandermonde = np.delete(deflection_vandermonde,3,0)
        print "7x7 CARTESIAN????"
        print cartesian_vandermonde

        cartesian_vandermonde_inv = self._left_inverse(cartesian_vandermonde)
        #cartesian_vandermonde_inv = inv(cartesian_vandermonde)
        transform_matrix = np.dot(cartesian_vandermonde_inv, deflection_vandermonde)

        print "TRANSFORM MATRIX!!!! SQUEEEEEE"
        print transform_matrix





        upper_distances = [distance for (deflection, distance) in upper_points.items()]
        upper_deflections = [deflection for (deflection, distance) in upper_points.items()]
        lower_distances = [distance for (deflection, distance) in lower_points.items()]
        lower_deflections = [deflection for (deflection, distance) in lower_points.items()]
        lower_height = 0


        #np.concatenate axis=1 is colum concatenate, axis=0 is concatenating a new row
        #Bring it all together into 8x3 matricies
        upper_3d_distances = np.concatenate((upper_distances, [(upper_height,)]*4), axis=1)
        lower_3d_distances = np.concatenate((lower_distances, [(lower_height,)]*4), axis=1)
        distances_3d = np.concatenate((upper_3d_distances, lower_3d_distances), axis=0)

        #concatenate the height and constant to finish the output side of the array
        upper_3d_deflections= np.concatenate((upper_deflections, [(upper_height,)]*4), axis=1)
        lower_3d_deflections= np.concatenate((lower_deflections, [(lower_height,)]*4), axis=1)
        deflections_3d = np.concatenate((upper_3d_deflections, lower_3d_deflections), axis=0)

        #Try offsets before and after, so transforms are around 0
        G_matrix = self._create_calibration_matricies(deflections_3d, distances_3d)
        self._transform = G_matrix
        return G_matrix

    def _get_xyz_centroids(self, input_matrix):
        '''Take set of points in [(dx,dy,z), (dx,dy,z), (...)]
            returns matrix [x_max_centroid, y_max_centroid, z_max_centroid, avg, x_min_centroid, y_min_centroid, z_min_centroid]
            where each centroid contains (x,y,z)'''
        centroids=[]
        x_cnt = 0
        y_cnt = 0
        z_cnt = 0

        for (i,xyz) in enumerate(input_matrix):
            x_cnt = x_cnt + xyz[0]
            y_cnt = y_cnt + xyz[1]
            z_cnt = z_cnt + xyz[2]

        rows = i+1
        x_avg = x_cnt/rows
        y_avg = y_cnt/rows
        z_avg = z_cnt/rows
        avg = [x_avg, y_avg, z_avg]
        x_max = [0,0,0]
        y_max = [0,0,0]
        z_max = [0,0,0]
        x_min = [0,0,0]
        y_min = [0,0,0]
        z_min = [0,0,0]

        for xyz in input_matrix:
            if xyz[0] > avg[0]:
                x_max[0] += (xyz[0] - x_avg)/4
                x_max[1] += (xyz[1] - y_avg)/4
                x_max[2] += (xyz[2] - z_avg)/4
            else:                  
                x_min[0] += (xyz[0] - x_avg)/4
                x_min[1] += (xyz[1] - y_avg)/4
                x_min[2] += (xyz[2] - z_avg)/4
            if xyz[1] > avg[1]:    
                y_max[0] += (xyz[0] - x_avg)/4
                y_max[1] += (xyz[1] - y_avg)/4
                y_max[2] += (xyz[2] - z_avg)/4
            else:                  
                y_min[0] += (xyz[0] - x_avg)/4
                y_min[1] += (xyz[1] - y_avg)/4
                y_min[2] += (xyz[2] - z_avg)/4
            if xyz[2] > avg[2]:    
                z_max[0] += (xyz[0] - x_avg)/4
                z_max[1] += (xyz[1] - y_avg)/4
                z_max[2] += (xyz[2] - z_avg)/4
            else:                  
                z_min[0] += (xyz[0] - x_avg)/4
                z_min[1] += (xyz[1] - y_avg)/4
                z_min[2] += (xyz[2] - z_avg)/4
        centroids=[x_max, x_min, y_max, y_min, z_max, z_min, avg]
        return centroids

    def _create_calibration_matricies(self, deflection_matrix, cartesian_matrix):

        #Get the centroids from the vertex points, and remove any offsets
        deflection_centroids = self._get_xyz_centroids(deflection_matrix)
        [xd_max, xd_min, yd_max, yd_min, zd_max, zd_min, d_avg] = deflection_centroids
        cartesian_centroids = self._get_xyz_centroids(cartesian_matrix)
        [x_max, x_min, y_max, y_min, z_max, z_min, avg] = cartesian_centroids

        x = 0
        y = 1
        z = 2

        # cheating a bit, the 0 offset is given in average of all centroids
        d_bx = d_avg[x]
        d_by = d_avg[y]
        d_bz = d_avg[z]
        bx = avg[x]
        by = avg[y]
        bz = avg[z]

        self._deflection_offsets = np.matrix([(d_bx, d_by, d_bz)])
        self._cartesian_offsets = np.matrix([(bx, by, bz)])

        #Create the Vandermonde matrix with the centroids:
        vandermonde_c_matrix = [(x_max[0], x_max[1], x_max[2], x_max[0]**2, x_max[1]**2, x_max[2]**2, 1),
                               (x_min[0], x_min[1], x_min[2], x_min[0]**2, x_min[1]**2, x_min[2]**2, 1),
                               (y_max[0], y_max[1], y_max[2], y_max[0]**2, y_max[1]**2, y_max[2]**2, 1),
                               (y_min[0], y_min[1], y_min[2], y_min[0]**2, y_min[1]**2, y_min[2]**2, 1),
                               (z_max[0], z_max[1], z_max[2], z_max[0]**2, z_max[1]**2, z_max[2]**2, 1),
                               (z_min[0], z_min[1], z_min[2], z_min[0]**2, z_min[1]**2, z_min[2]**2, 1),
                               (0, 0, 0, 0, 0, 0, 1)]

        vandermonde_d_matrix = [(xd_max[0], xd_max[1], xd_max[2], xd_max[0]**2, xd_max[1]**2, xd_max[2]**2, 1),
                               (xd_min[0], xd_min[1], xd_min[2], xd_min[0]**2, xd_min[1]**2, xd_min[2]**2, 1),
                               (yd_max[0], yd_max[1], yd_max[2], yd_max[0]**2, yd_max[1]**2, yd_max[2]**2, 1),
                               (yd_min[0], yd_min[1], yd_min[2], yd_min[0]**2, yd_min[1]**2, yd_min[2]**2, 1),
                               (zd_max[0], zd_max[1], zd_max[2], zd_max[0]**2, zd_max[1]**2, zd_max[2]**2, 1),
                               (zd_min[0], zd_min[1], zd_min[2], zd_min[0]**2, zd_min[1]**2, zd_min[2]**2, 1),
                               (0, 0, 0, 0, 0, 0, 1)]

        c_matrix_inv = self._left_inverse(vandermonde_c_matrix)
        vandermonde_t_matrix = np.dot(c_matrix_inv, vandermonde_d_matrix)

        #print "Vander C Matrix:"
        #print np.matrix(vandermonde_c_matrix)
        #print "Vander D Matrix:"
        #print np.matrix(vandermonde_d_matrix)
        print "VANDERMONDE MATRIX:"
        print vandermonde_t_matrix

        return vandermonde_t_matrix


    def _left_inverse(self,matrix):
        ''' leftInverse(A) = inverse(A_t*A)*A_t
            This should take the moore-penrose inverse'''

        transposed_matrix = np.transpose(matrix)
        sq_helper = np.dot(transposed_matrix,matrix)
        inverse_helper = np.linalg.inv(sq_helper)
        left_inverse = np.dot(inverse_helper,transposed_matrix)
        return left_inverse

    def transform(self,xyz_cartesian):
        '''THE transform function, provide (x,y,z) and you'll get (xd, yd, z)'''
        xyz_point = (np.matrix([xyz_cartesian]) - self._cartesian_offsets).tolist()
        xyz_v=[(xyz_point[0][0], xyz_point[0][1], xyz_point[0][2], xyz_point[0][0]**2, xyz_point[0][1]**2, xyz_point[0][2]**2, 1)]
        deflections = np.dot(xyz_v,self._transform)
        deflections = [(deflections[0][0], deflections[0][1], deflections[0][2])] + self._deflection_offsets

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

        self._offset_params = {'mx':0,'my':0,'bx':0,'by':0}
        (upper_deflections, upper_cartesian) = zip(*upper_points)
        (lower_deflections, lower_cartesian) = zip(*lower_points)

        centroids = self._get_centroids(upper_deflections, lower_deflections)
        self._create_offset_functions(centroids, upper_height)

        upper_deflections = self._remove_xy_offsets(upper_deflections, upper_height)
        lower_deflections = self._remove_xy_offsets(lower_deflections, 0)
        upper_points = zip(upper_deflections, upper_cartesian)
        lower_points = zip(lower_deflections, lower_cartesian)

        #take [opposite points][Deflections][x or y] as each of the options
        deflection_scale = (
            (upper_points[0][0][0] - upper_points[2][0][0]) / (lower_points[0][0][0] - lower_points[2][0][0]),
            (upper_points[0][0][1] - upper_points[2][0][1]) / (lower_points[0][0][1] - lower_points[2][0][1])
            )

        self._lower_points = lower_points
        self._upper_points = [(self._scale_point(deflection, deflection_scale), distance) for (deflection, distance) in lower_points]

        self._get_transforms()
        self._cache = {}

    def _remove_xy_offsets(self,xy_list, z_height):
        '''Returns (x,y,z) with the centroid offset removed to bring deflections between -0.5 to 0.5'''
        result_matrix = []
        for xy in xy_list:
            x_offset = z_height*self._offset_params['mx'] + self._offset_params['bx']
            y_offset = z_height*self._offset_params['my'] + self._offset_params['by']
            x_prime = xy[0]-x_offset
            y_prime = xy[1]-y_offset
            result_matrix.append((x_prime,y_prime))
        return result_matrix

    def _add_xy_offset(self,xyz):
        '''Returns (x,y,z) with the centroid offset added to bring deflections between 0-1'''
        z_height = xyz[2]
        x_offset = z_height*self._offset_params['mx'] + self._offset_params['bx']
        y_offset = z_height*self._offset_params['my'] + self._offset_params['by']
        x_prime = xyz[0]+x_offset
        y_prime = xyz[1]+y_offset
        return (x_prime,y_prime)

    def _get_centroids(self, upper_deflections, lower_deflections):
        '''Take set of points in [((x,y)...(x,y)), ... ] and return [(x,y),(...)] centroid list
            - This may also just be the det(points_matrix) '''
        points = [upper_deflections,lower_deflections]
        centroids=[]
        for pointset in points:
            x = 0
            y = 0
            for (i, xyz) in enumerate(pointset):
                x=x+xyz[0]
                y=y+xyz[1]
            total=i+1
            centroids.append((x/total,y/total))
        return centroids

    def _create_offset_functions(self, centroids, height):
        '''Takes the matrix of (upper;lower) centroids and creates the offset transform
           to satisfy the offset part of calibration
           
           This creates the ol y=mx+b but it derives dx=m*z+b'''

        top_centroid = centroids[0]
        bottom_centroid = centroids[1]
        top_x = top_centroid[0]
        top_y = top_centroid[1]
        top_z = height
        bottom_x = bottom_centroid[0]
        bottom_y = bottom_centroid[1]
        bottom_z = 0

        mx = (top_x - bottom_x)/(top_z - bottom_z)
        my = (top_y - bottom_y)/(top_z - bottom_z)
        bx = top_x - bottom_x*mx
        by = top_y - bottom_y*my

        self._offset_params['mx']=mx
        self._offset_params['my']=my
        self._offset_params['bx']=bx
        self._offset_params['by']=by

        return

    def _scale_point(self, point, scale):
        x, y = point
        sx = x * scale[0]
        sy = y * scale[1]
        return (sx, sy)

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
        (x2, y2) = self._add_xy_offset((x1, y1, z))
        if x2 >= 0.0 and x2 <= 1.0 and y2 >= 0.0 and y2 <= 1.0:
            return (x2, y2)
        else:
            logger.warning("Bounds of printer exceeded: %s,%s" % (x, y))
            adjusted_x = min(1.0, max(0.0, x1))
            adjusted_y = min(1.0, max(0.0, y1))
            #re-add offsets
            return(adjusted_translated_x, adjusted_translated_y)

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

    #print "Upper points:"
    #print upper_points

    #print "Lower Points:"
    #print lower_points
    #
    example_xyz = (0.0,0.0,0.0)

    scale = 1
    #HomoTransformer = HomogenousTransformer(scale, height, lower_points_real, upper_points_real)
    if (True):
        print "LinTransformerMade"
        LinTransformer=LinearAlgebraTransformer(height ,lower_points_real ,upper_points_real)
        #LinTransformer=LinearAlgebraTransformer(height ,lower_points, upper_points)
        print "Input to transform function: {0}".format(example_xyz)
        deflections = LinTransformer.transform(example_xyz)
        print "Deflections at 0mm centered: {0}".format(deflections)

        example_xyz = (20.0,20.0,10.0)
        print example_xyz
        deflections = LinTransformer.transform(example_xyz)
        print "Deflections for {0} = {1}".format(example_xyz, deflections)
        example_xyz = (-20.0,-20.0,10.0)
        deflections = LinTransformer.transform(example_xyz)
        print "Deflections for {0} = {1}".format(example_xyz, deflections)
        example_xyz = (-20.0,-20.0,50.0)
        deflections = LinTransformer.transform(example_xyz)
        print "Deflections for {0} = {1}".format(example_xyz, deflections)

        example_xyz = (10.0,10.0,20.0)
        deflections = LinTransformer.transform(example_xyz)
        print "Deflections for {0} = {1}".format(example_xyz, deflections)
        example_xyz = (10.0,10.0,40.0)
        deflections = LinTransformer.transform(example_xyz)
        print "Deflections for {0} = {1}".format(example_xyz, deflections)
        example_xyz = (10.0,10.0,80.0)
        deflections = LinTransformer.transform(example_xyz)
        print "Deflections for {0} = {1}".format(example_xyz, deflections)

        example_xyz = (0.0,0.0,50.0)
        deflections = LinTransformer.transform(example_xyz)
        print "Deflections after 50mm centered: {0}".format(deflections)

        import time
        now=time.time()
        iterations=100000
        for i in range(1,iterations):
            pass
            deflections = LinTransformer.transform(example_xyz)
        after=time.time()
        print "{0} iterations of transform took {1} seconds".format(iterations,after-now)

    
