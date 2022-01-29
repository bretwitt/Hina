cmake_minimum_required(VERSION 3.0.0)
project(Hina VERSION 0.1.0)

add_library (Hina STATIC
                src/core/solver/WheelForceSolver.h  src/core/solver/WheelForceSolver.cpp
                src/core/environment/soil/SoilPatch.h src/core/environment/soil/SoilPatch.cpp                                   
                src/core/environment/terrain/TerrainPatch.h
                src/core/solver/BekkerForceSolver/BekkerForceSolver.h src/core/solver/BekkerForceSolver/BekkerForceSolver.cpp 
                src/core/solver/WheelForceSolver.h
                src/core/body/Body.h src/core/body/Wheel.h
            )

add_executable(Hina_Driver Hina.cpp)
target_link_libraries(Hina_Driver Hina gsl)