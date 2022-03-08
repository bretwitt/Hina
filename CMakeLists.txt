cmake_minimum_required(VERSION 3.0.0)
project(Hina)


add_library (Hina STATIC
                src/core/solver/WheelForceSolver.h  src/core/solver/WheelForceSolver.cpp
                src/core/environment/soil/SoilPatch.h src/core/environment/soil/SoilPatch.cpp                                   
                src/core/environment/terrain/TerrainPatch.h
                src/core/solver/BekkerForceSolver/BekkerForceSolver.h src/core/solver/BekkerForceSolver/BekkerForceSolver.cpp 
                src/core/solver/WheelForceSolver.h
                src/core/body/Body.h src/core/body/Wheel.h
            )

        
include_directories(/usr/include/python3.8)

add_executable(Hina_Driver Hina.cpp)

target_compile_features(Hina_Driver PRIVATE cxx_std_17)

target_link_libraries(Hina_Driver Hina gsl gslcblas)