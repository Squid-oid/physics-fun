# physics-fun
A little project to build and visualize a physics engine

# Todo

Phase 1 (Shapes in Box);
1. Create Rectangle W_object and test it's collision handling
2. a) Create a tri object
|. b) Create a Generic Conv_poly object composed of multiple tris
4. Replace Rectangle W_object implementation with Conv_poly object backend

Phase 2 (Good Shapes in Box);
1. Implement multiple types of object materials which give different collsion properties (elastic, semi elastic, force killing different weights etc.)
2. Create Test Suite for all previously created work
3. Implement External Forces, starting with gravity, maybe also Internal Value, Coordinate and Velocity dependant forces ie Electromagnetic
4. Make one of the Shapes Controllable to provide 'player' interaction

Phase 3 (3d Spheres);
[Create 3d spheres in box simulation using what I've Learned]

# Completed
Phase 0 (Balls in Box); COMPLETE
1. ~~Add proper comments to code already present~~
2. ~~Reevaluate the necessity of the CVector Vector handling~~ [Removed CVector usage]
3. ~~Change refresh() actions over to using the TimeFuncs.timestep() function, as well as being tied to the system clock instead of framerate~~
4. ~~Implement and test bsic ball on ball collision~~
5. ~~Allow free selection of Ball Size~~
6. ~~[Maybe] implement refresh() new collision rechecking to allow a Ball object to continue bouncing multiple times in one refresh~~, [Maybe] ~~implement collision based physics ticks rather than frame based ones, ie extrapolate object paths until first collision/(lazily possible collision using event horizon) and do a physics refresh exactly at that point to reduce tick frequency.~~ Second part Unfeasible at this stage
7. ~~Fix Basic Dynamic Sprite Generation for Barriers and Balls~~
8. ~~Fix physics bug (Likely in collision)~~        RESOLVED - Double counting of certain offsets was issue
9. ~~Work on gfx bugs~~
10. ~~Different weights~~
