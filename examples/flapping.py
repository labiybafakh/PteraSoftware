import pterasoftware as ps
import numpy as np


ornibibot = ps.geometry.airplane.Airplane(
    name="Flapping Bird",
    Cg_GP1_CgP1=(0.0, 0.0, 0.0),
    weight=0.38,
    s_ref=None,
    c_ref=None,
    b_ref=None,
    wings=[
        ps.geometry.wing.Wing(
            name="Main Wing",
            Ler_Gs_Cgs=(0.0, 0.03, 0.0),
            angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
            symmetric=True,
            mirror_only=False,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            num_chordwise_panels=6,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=8,
                    chord=0.32,
                    Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.001,
                    control_surface_deflection=0.0,
                    spanwise_spacing="cosine",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="AS6095",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=8,
                    chord=0.32,
                    Lp_Wcsp_Lpp=(0.0, 0.15, 0.0),
                    angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.001,
                    control_surface_deflection=0.0,
                    spanwise_spacing="cosine",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="AS6095",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=8,
                    chord=0.32,
                    Lp_Wcsp_Lpp=(0.0, 0.3, 0.0),
                    angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.001,
                    control_surface_deflection=0.0,
                    spanwise_spacing="cosine",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="AS6095",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=None,
                    chord=0.001,
                    Lp_Wcsp_Lpp=(0.0, 0.75, 0.0),
                    angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.001,
                    control_surface_deflection=0.0,
                    spanwise_spacing=None,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="AS6095",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
            ],
        ),
        ps.geometry.wing.Wing(
            name="Tail",
            Ler_Gs_Cgs=(0.5, 0.0, 0.0),
            angles_Gs_to_Wn_ixyz=(0.0, -5.0, 0.0),
            symmetric=True,
            mirror_only=False,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            num_chordwise_panels=6,
            chordwise_spacing="uniform",
            # Define this wing's root wing cross section.
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=8,
                    chord=0.3,
                    Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.001,
                    control_surface_deflection=0.0,
                    spanwise_spacing="cosine",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="2032c",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=None,
                    chord=0.001,
                    Lp_Wcsp_Lpp=(0.25, 0.2, 0.0),
                    angles_Wcsp_to_Wcs_ixyz=(0.0, -20.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.001,
                    control_surface_deflection=0.0,
                    spanwise_spacing=None,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="2032c",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
            ],
        ),
    ],
)

print(f"   Created airplane: {ornibibot.name}")
print(f"   Number of wings: {len(ornibibot.wings)}")
print(f"   Wing span: {ornibibot.wings[0].span:.2f} m")



# Now define the main wing's WingCrossSections' WingCrossSectionMovements for all 4 sections.
main_wing_root_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=ornibibot.wings[0].wing_cross_sections[0],
        ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )
)
main_wing_mid1_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=ornibibot.wings[0].wing_cross_sections[1],
        ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(60.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(0.33, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )
)
main_wing_mid2_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=ornibibot.wings[0].wing_cross_sections[2],
        ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )
)
main_wing_tip_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=ornibibot.wings[0].wing_cross_sections[3],
        ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )
)

# The main Wing was defined to have symmetric=True, mirror_only=False, and with a
# symmetry plane offset non coincident with the Wing's axes yz plane. Therefore,
# that Wing had type 5 symmetry (see the Wing class documentation for more details on
# symmetry types). Therefore, it was actually split into two Wings, the with the
# second Wing being a reflected version of the first. Therefore, we need to define a
# WingMovement for this reflected Wing. To start, we'll first define the reflected
# main wing's WingCrossSections' WingCrossSectionMovements for all sections.
reflected_main_wing_root_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=ornibibot.wings[1].wing_cross_sections[0],
        ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )
)
reflected_main_wing_mid1_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=ornibibot.wings[1].wing_cross_sections[1],
        ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )
)
reflected_main_wing_mid2_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=ornibibot.wings[1].wing_cross_sections[2],
        ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )
)
reflected_main_wing_tip_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=ornibibot.wings[1].wing_cross_sections[3],
        ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )
)

# Now define the v tail's root and tip WingCrossSections' WingCrossSectionMovements.
v_tail_root_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=ornibibot.wings[2].wing_cross_sections[0],
        ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )
)
v_tail_tip_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=ornibibot.wings[2].wing_cross_sections[1],
        ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )
)

# Now define the main wing's WingMovement, the reflected main wing's WingMovement and
# the v tail's WingMovement.
main_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=ornibibot.wings[0],
    wing_cross_section_movements=[
        main_wing_root_wing_cross_section_movement,
        main_wing_mid1_wing_cross_section_movement,
        main_wing_mid2_wing_cross_section_movement,
        main_wing_tip_wing_cross_section_movement,
    ],
    ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
    periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
    spacingLer_Gs_Cgs=("sine", "sine", "sine"),
    phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
    ampAngles_Gs_to_Wn_ixyz=(60.0, 0.0, 0.0),  # Reduced amplitude for smoother motion
    periodAngles_Gs_to_Wn_ixyz=(0.33, 0.0, 0.0),   # Synchronized with reflected wing
    spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
    phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
)
reflected_main_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=ornibibot.wings[1],
    wing_cross_section_movements=[
        reflected_main_wing_root_wing_cross_section_movement,
        reflected_main_wing_mid1_wing_cross_section_movement,
        reflected_main_wing_mid2_wing_cross_section_movement,
        reflected_main_wing_tip_wing_cross_section_movement,
    ],
    ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
    periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
    spacingLer_Gs_Cgs=("sine", "sine", "sine"),
    phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
    ampAngles_Gs_to_Wn_ixyz=(60.0, 0.0, 0.0),  # Synchronized amplitude
    periodAngles_Gs_to_Wn_ixyz=(0.33, 0.0, 0.0),  # Same period
    spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
    phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
)
v_tail_movement = ps.movements.wing_movement.WingMovement(
    base_wing=ornibibot.wings[2],
    wing_cross_section_movements=[
        v_tail_root_wing_cross_section_movement,
        v_tail_tip_wing_cross_section_movement,
    ],
    ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
    periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
    spacingLer_Gs_Cgs=("sine", "sine", "sine"),
    phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
    ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    periodAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
    phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
)

# Delete the extraneous pointers to the WingCrossSectionMovements, as these are now
# contained within the WingMovements. This is optional, but it can make debugging
# easier.
del main_wing_root_wing_cross_section_movement
del main_wing_mid1_wing_cross_section_movement
del main_wing_mid2_wing_cross_section_movement
del main_wing_tip_wing_cross_section_movement
del reflected_main_wing_root_wing_cross_section_movement
del reflected_main_wing_mid1_wing_cross_section_movement
del reflected_main_wing_mid2_wing_cross_section_movement
del reflected_main_wing_tip_wing_cross_section_movement
del v_tail_root_wing_cross_section_movement
del v_tail_tip_wing_cross_section_movement

# Now define the example airplane's AirplaneMovement.
airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=ornibibot,
    wing_movements=[main_wing_movement, reflected_main_wing_movement, v_tail_movement],
    ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
    periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
    spacingCg_GP1_CgP1=("sine", "sine", "sine"),
    phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
)

# Delete the extraneous pointers to the WingMovements.
del main_wing_movement
del reflected_main_wing_movement
del v_tail_movement

# Define a new OperatingPoint.
example_operating_point = ps.operating_point.OperatingPoint(
    rho=1.225, vCg__E=4.2, alpha=1.0, beta=0.0, externalFX_W=0.0, nu=15.06e-6
)

# Define the operating point's OperatingPointMovement.
operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
    base_operating_point=example_operating_point, periodVCg__E=0.0, spacingVCg__E="sine"
)

# Delete the extraneous pointer.
del example_operating_point

# Define the Movement. This contains the AirplaneMovement and the
# OperatingPointMovement.
movement = ps.movements.movement.Movement(
    airplane_movements=[airplane_movement],
    operating_point_movement=operating_point_movement,
    delta_time=0.001,
    num_cycles=None,
    num_chords=None,
    num_steps=100,
)

# Delete the extraneous pointers.
del airplane_movement
del operating_point_movement

# Define the UnsteadyProblem.
example_problem = ps.problems.UnsteadyProblem(
    movement=movement,
)

# Define a new solver. The available solver classes are
# SteadyHorseshoeVortexLatticeMethodSolver, SteadyRingVortexLatticeMethodSolver,
# and UnsteadyRingVortexLatticeMethodSolver. We'll create an
# UnsteadyRingVortexLatticeMethodSolver, which requires a UnsteadyProblem.
example_solver = (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=example_problem,
    )
)

# Delete the extraneous pointer.
del example_problem

# Run the solver.
example_solver.run(
    prescribed_wake=True,
    show_progress=True,
)

# Call the animate function on the solver. This produces a GIF of the wake being
# shed. The GIF is saved in the same directory as this script. Press "q",
# after orienting the view, to begin the animation.
ps.output.animate(
    unsteady_solver=example_solver,
    scalar_type="lift",
    show_wake_vortices=True,
    save=False,
)
