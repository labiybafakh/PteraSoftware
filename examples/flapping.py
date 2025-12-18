import pterasoftware as ps
import numpy as np


airplane = ps.geometry.Airplane(
    name="Flapping Bird",
    wings=[
        ps.geometry.Wing(
            name="Main Wing",
            symmetric=True,
            num_chordwise_panels=8,
            chordwise_spacing="cosine",
            wing_cross_sections=[
                # Root cross section
                ps.geometry.WingCrossSection(
                    x_le=0.0,  # Leading edge x position
                    y_le=0.0,  # Leading edge y position (root)
                    z_le=0.0,  # Leading edge z position
                    chord=0.32,  # Chord length in meters
                    twist=0.0,  # Twist angle in degrees
                    airfoil=ps.geometry.Airfoil(name="as6095"),
                    num_spanwise_panels=4,
                ),
                # Mid cross section
                ps.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=0.15,  
                    z_le=0.0,
                    chord=0.32,
                    twist=0.0,
                    airfoil=ps.geometry.Airfoil(name="flat_plate"),
                    num_spanwise_panels=4,
                ),
                # Tip cross section
                ps.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=0.75,  
                    z_le=0.0,
                    chord=0.01,
                    twist=0.0,
                    airfoil=ps.geometry.Airfoil(name="flat_plate"),
                    num_spanwise_panels=4,
                ),
            ],
        ),
        ps.geometry.Wing(
            name="V-Tail",
            x_le=6.75,
            z_le=0.25,
            num_chordwise_panels=6,
            chordwise_spacing="uniform",
            symmetric=True,
            # Define this wing's root wing cross section.
            wing_cross_sections=[
                ps.geometry.WingCrossSection(
                    chord=0.01,
                    # Give the root wing cross section an airfoil.
                    airfoil=ps.geometry.Airfoil(
                        name="flat_plate",
                    ),
                    twist=5.0,  # Give the root wing cross section a twist.
                ),
                # Define the wing's tip wing cross section.
                ps.geometry.WingCrossSection(
                    x_le=0.5,
                    y_le=2.0,
                    z_le=1.0,
                    chord=1.0,
                    twist=-5.0,  # Give the tip wing cross section an airfoil.
                    airfoil=ps.geometry.Airfoil(
                        name="flat_plate",
                    ),
                ),
            ],
        ),
    ],
)

print(f"   Created airplane: {airplane.name}")
print(f"   Number of wings: {len(airplane.wings)}")
print(f"   Wing span: {airplane.wings[0].span:.2f} m")



# Now define the main wing's root wing cross section's movement. Cross sections can
# move in three ways: sweeping, pitching, and heaving. Sweeping is defined as the
# relative rotation of this wing cross section's leading edge to its preceding wing
# cross section's leading edge about the airplane's body x axis. Pitching is defined
# as the relative rotation of this wing cross section's leading edge to the preceding
# wing cross section's leading edge about the body y axis. Heaving is defined as the
# relative rotation of this wing cross section's leading edge to the preceding wing
# cross section's leading edge about the body z axis. The sign of all rotations is
# determined via the right-hand-rule.
main_wing_root_wing_cross_section_movement = ps.movement.WingCrossSectionMovement(
    # Provide the base cross section.
    base_wing_cross_section=airplane.wings[0].wing_cross_sections[0],
    # Define the sweeping amplitude. This value is in degrees. As this is the first
    # wing cross section, this must be 0.0 degrees, which is the default value.
    sweeping_amplitude=0.0,
    # Define the sweeping period. This value is in seconds. As this is the first wing
    # cross section, this must be 0.0 seconds, which is the default value.
    sweeping_period=0.0,
    # Define the time step spacing of the sweeping. This is "sine" by default. The
    # options are "sine" and "uniform".
    sweeping_spacing="sine",
    # Define the pitching amplitude. This value is in degrees. As this is the first
    # wing cross section, this must be 0.0 degrees, which is the default value.
    pitching_amplitude=0.0,
    # Define the pitching period. This value is in seconds. As this is the first wing
    # cross section, this must be 0.0 seconds, which is the default value.
    pitching_period=0.0,
    # Define the time step spacing of the pitching. This is "sine" by default. The
    # options are "sine" and "uniform".
    pitching_spacing="sine",
    # Define the heaving amplitude. This value is in degrees. As this is the first
    # wing cross section, this must be 0.0 degrees, which is the default value.
    heaving_amplitude=0.0,
    # Define the heaving period. This value is in seconds. As this is the first wing
    # cross section, this must be 0.0 seconds, which is the default value.
    heaving_period=0.0,
    # Define the time step spacing of the heaving. This is "sine" by default. The
    # options are "sine" and "uniform".
    heaving_spacing="sine",
)

# Define the main wing's tip wing cross section's movement.
main_wing_tip_wing_cross_section_movement = ps.movement.WingCrossSectionMovement(
    base_wing_cross_section=airplane.wings[0].wing_cross_sections[2],
    sweeping_amplitude=60.0,
    sweeping_period=0.33,
    sweeping_spacing="sine",
    pitching_amplitude=15.0,
    pitching_period=1.0,
    pitching_spacing="sine",
    heaving_amplitude=0.0,
    heaving_period=0.0,
    heaving_spacing="sine",
)

# Define the v-tail's root wing cross section's movement. This wing will be static,
# so the movement attributes can be excluded, and the default values will suffice.
v_tail_root_wing_cross_section_movement = ps.movement.WingCrossSectionMovement(
    base_wing_cross_section=airplane.wings[1].wing_cross_sections[0],
)

# Define the v-tail's root wing cross section's movement. This wing will be static,
# so the movement attributes can be excluded, and the default values will suffice.
v_tail_tip_wing_cross_section_movement = ps.movement.WingCrossSectionMovement(
    base_wing_cross_section=airplane.wings[1].wing_cross_sections[1],
)

# Now define the main wing's movement. In addition to their wing cross sections'
# relative movements, wings' leading edge positions can move as well.
main_wing_movement = ps.movement.WingMovement(  # Define the base wing object.
    base_wing=airplane.wings[0],
    # Add the list of wing cross section movement objects.
    wing_cross_sections_movements=[
        main_wing_root_wing_cross_section_movement,
        main_wing_tip_wing_cross_section_movement,
    ],
    # Define the amplitude of the leading edge's change in x position. This value is
    # in meters. This is set to 0.0 meters, which is the default value.
    x_le_amplitude=0.0,
    # Define the period of the leading edge's change in x position. This is set to
    # 0.0 seconds, which is the default value.
    x_le_period=0.0,
    # Define the time step spacing of the leading edge's change in x position. This
    # is "sine" by default. The options are "sine" and "uniform".
    x_le_spacing="sine",
    # Define the amplitude of the leading edge's change in y position. This value is
    # in meters. This is set to 0.0 meters, which is the default value.
    y_le_amplitude=0.0,
    # Define the period of the leading edge's change in y position. This is set to
    # 0.0 seconds, which is the default value.
    y_le_period=0.0,
    # Define the time step spacing of the leading edge's change in y position. This
    # is "sine" by default. The options are "sine" and "uniform".
    y_le_spacing="sine",
    # Define the amplitude of the leading edge's change in z position. This value is
    # in meters. This is set to 0.0 meters, which is the default value.
    z_le_amplitude=0.0,
    # Define the period of the leading edge's change in z position. This is set to
    # 0.0 seconds, which is the default value.
    z_le_period=0.0,
    # Define the time step spacing of the leading edge's change in z position. This
    # is "sine" by default. The options are "sine" and "uniform".
    z_le_spacing="sine",
)

# Delete the extraneous wing cross section movement objects, as these are now
# contained within the wing movement object. This is unnecessary, but it can make
# debugging easier.
del main_wing_root_wing_cross_section_movement
del main_wing_tip_wing_cross_section_movement

# Make the v-tail's wing movement object.
v_tail_movement = ps.movement.WingMovement(  # Define the base wing object.
    base_wing=airplane.wings[1],
    # Add the list of wing cross section movement objects.
    wing_cross_sections_movements=[
        v_tail_root_wing_cross_section_movement,
        v_tail_tip_wing_cross_section_movement,
    ],
)

# Delete the extraneous wing cross section movement objects, as these are now
# contained within the wing movement object. This is unnecessary, but it can make
# debugging easier.
del v_tail_root_wing_cross_section_movement
del v_tail_tip_wing_cross_section_movement

# Now define the airplane's movement object. In addition to their wing's and wing
# cross sections' relative movements, airplane's reference positions can move as well.
airplane_movement = ps.movement.AirplaneMovement(  # Define the base airplane object.
    base_airplane=airplane,  # Add the list of wing movement objects.
    wing_movements=[main_wing_movement, v_tail_movement],
    # Define the amplitude of the reference position's change in x position. This
    # value is in meters. This is set to 0.0 meters, which is the default value.
    x_ref_amplitude=0.0,
    # Define the period of the reference position's change in x position. This value
    # is in seconds. This is set to 0.0 seconds, which is the default value.
    x_ref_period=0.0,
    # Define the time step spacing of the reference position's change in x position.
    # This is "sine" by default. The options are "sine" and "uniform".
    x_ref_spacing="sine",
    # Define the amplitude of the reference position's change in y position. This
    # value is in meters. This is set to 0.0 meters, which is the default value.
    y_ref_amplitude=0.0,
    # Define the period of the reference position's change in y position. This value
    # is in seconds. This is set to 0.0 seconds, which is the default value.
    y_ref_period=0.0,
    # Define the time step spacing of the reference position's change in y position.
    # This is "sine" by default. The options are "sine" and "uniform".
    y_ref_spacing="sine",
    # Define the amplitude of the reference position's change in z position. This
    # value is in meters. This is set to 0.0 meters, which is the default value.
    z_ref_amplitude=0.0,
    # Define the period of the reference position's change in z position. This value
    # is in seconds. This is set to 0.0 seconds, which is the default value.
    z_ref_period=0.0,
    # Define the time step spacing of the reference position's change in z position.
    # This is "sine" by default. The options are "sine" and "uniform".
    z_ref_spacing="sine",
)

# Delete the extraneous wing movement objects, as these are now contained within the
# airplane movement object.
del main_wing_movement
del v_tail_movement

# Define a new operating point object. This defines the state at which the airplane
# object is operating.
example_operating_point = ps.operating_point.OperatingPoint(
    # Define the density of the fluid the airplane is flying in. This defaults to
    # 1.225 kilograms per meters cubed.
    density=1.225,
    # Define the angle of sideslip the airplane is experiencing. This defaults to 0.0
    # degrees.
    beta=0.0,
    # Define the freestream velocity at which the airplane is flying. This defaults
    # to 10.0 meters per second.
    velocity=10.0,
    # Define the angle of attack the airplane is experiencing. This defaults to 5.0
    # degrees.
    alpha=1.0,
    # Define the kinematic viscosity of the air in meters squared per second. This
    # defaults to 15.06e-6 meters squared per second, which corresponds to an air
    # temperature of 20 degrees Celsius.
    nu=15.06e-6,
)

# Define the operating point's movement. The operating point's velocity can change
# with respect to time.
operating_point_movement = ps.movement.OperatingPointMovement(
    # Define the base operating point object.
    base_operating_point=example_operating_point,
    # Define the amplitude of the velocity's change in time. This value is set to 0.0
    # meters per second, which is the default value.
    velocity_amplitude=0.0,
    # Define the period of the velocity's change in time. This value is set to 0.0
    # seconds, which is the default value.
    velocity_period=0.0,
    # Define the time step spacing of the velocity's change in time. This is "sine"
    # by default. The options are "sine" and "uniform".
    velocity_spacing="sine",
)

# Define the movement object. This contains the airplane movement and the operating
# point movement.
movement = ps.movement.Movement(  # Add the airplane movement.
    airplane_movements=[airplane_movement],  # Add the operating point movement.
    operating_point_movement=operating_point_movement,
    # Leave the number of time steps and the length of each time step unspecified.
    # The solver will automatically set the length of the time steps so that the wake
    # ring vortices and the bound ring vortices have approximately the same area. The
    # solver will also determine if the geometry is static or not. If it is static,
    # the number of steps will be set such that the wake extends ten chord lengths
    # back from the main wing. If the geometry isn't static, the number of steps will
    # be set such that three periods of the slowest movement oscillation complete.
    num_steps=None,
    delta_time=None,
)

# Delete the extraneous airplane and operating point movement objects, as these are
# now contained within the movement object.
del airplane_movement
del operating_point_movement

# Define the unsteady example problem.
example_problem = ps.problems.UnsteadyProblem(
    movement=movement,
)

# Define a new solver. The available solver objects are the steady horseshoe vortex
# lattice method solver, the steady ring vortex lattice method solver, and the
# unsteady ring vortex lattice method solver.
example_solver = ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
    # Solvers just take in one attribute: the problem they are going to solve.
    unsteady_problem=example_problem,
)

# Delete the extraneous pointer to the problem as it is now contained within the
# solver.
del example_problem

# Run the example solver.
example_solver.run(
    # This parameter determines the detail of information that the solver's logger
    # will output while running. The options are, in order of detail and severity,
    # "Debug", "Info", "Warning", "Error", "Critical". The default value is "Warning".
    logging_level="Warning",
    # Use a prescribed wake model. This is faster, but may be slightly less accurate.
    prescribed_wake=True,
)

# Call the software's animate function on the solver. This produces a GIF of the wake
# being shed. The GIF is saved in the same directory as this script. Press "q",
# after orienting the view, to begin the animation.
ps.output.animate(  # Set the unsteady solver to the one we just ran.
    unsteady_solver=example_solver,
    # Tell the animate function to color the aircraft's wing panels with the local
    # lift coefficient. The valid arguments for this parameter are None, "induced drag",
    # "side force", or "lift".
    scalar_type="lift",
    # Tell the animate function to show the wake vortices. This value defaults to
    # False.
    show_wake_vortices=True,
    # Tell the animate function to not save the animation as file. This way,
    # the animation will still be displayed but not saved. This value defaults to
    # False.
    save=False,
)

# Compare the output you see with the expected outputs saved in the "docs/examples
# expected output" directory.
