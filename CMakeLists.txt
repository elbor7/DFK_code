cmake_minimum_required(VERSION 3.24)
project(PCC_Processing_Design_tool)

set(CMAKE_CXX_STANDARD 20)

add_executable( PCC_Processing_Design_tool
        # Management module
        src/main.cpp
        # Various support functions
        src/DCC_SupportFunctions.h
        # C++  class for PCC and related chains of special elements
        src/DCC_Objects/DCC_Objects.h
        # Processing module
        src/DCC_Processing/DCC_Processing.h
        src/DCC_Processing/Processing_Functions.h
        # Kinetic module
        src/DCC_Kinetic/DCC_Kinetic.h
        src/DCC_Kinetic/Kinetic_Functions.h
        # Characterisation module
        src/DCC_Characterisation/StructureCharacterisation.h
        src/DCC_Characterisation/TJsLab.h
        src/DCC_Characterisation/QPsLab.h
        src/DCC_Characterisation/DislocationsLab.h
        # Writer module
        src/DCC_Writer/DCC_Writer.h
        src/DCC_Writer/Writer_Functions.h
        # Subcomplex module
        src/DCC_Section/Planecut_Functions.h
        # Multiphysics module
        src/DCC_Multiphysics/DCC_Multiphysics.h
        src/DCC_Multiphysics/Multiphysics_Functions.h
        # Simulation Tasks here
        #src/tasks/task_macrocrack.cpp
        )
