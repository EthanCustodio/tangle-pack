import tanglepack
from tanglepack.examples import (
    HENON_K10,
    henon_map as _henon_map_factory,
    henon_map_inverse as _henon_map_inverse_factory,
)
import matplotlib.pyplot as plt


henon_map = _henon_map_factory(*HENON_K10)
henon_map_inverse = _henon_map_inverse_factory(*HENON_K10)


henon = tanglepack.DynamicalSystem(henon_map, henon_map_inverse)

initial_guess = [4, -4]

fp_solver = tanglepack.FixedPointSolver(henon)

fixed_point = fp_solver.construct_fixed_point(initial_guess)

print(f'The fixed point is: {fixed_point.coordinates[0]}')
print(f'The fixed point is type: {type(fixed_point)}')

man_maker = tanglepack.ManifoldInitializer(henon)

initial_segment = man_maker.get_initial_fundamental_segment(fixed_point, 0, 0, 'unstable')
initial_points = initial_segment.get_point_array()
print(f'Initial 3 Points: {initial_points}')


initial_segment.plot(marker='o', ms=10)

