##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=HW11_Tree
ConfigurationName      :=Debug
WorkspacePath          :=/home/zecophy/Desktop/numercialMethod/9821-Fall-2018-Team-5/HW11
ProjectPath            :=/home/zecophy/Desktop/numercialMethod/9821-Fall-2018-Team-5/HW11/HW11_Tree
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Chenyu Zhao
Date                   :=13/12/18
CodeLitePath           :=/home/zecophy/.codelite
LinkerName             :=/usr/bin/g++
SharedObjectLinkerName :=/usr/bin/g++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.i
DebugSwitch            :=-g 
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=$(IntermediateDirectory)/$(ProjectName)
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="HW11_Tree.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/ar rcu
CXX      := /usr/bin/g++
CC       := /usr/bin/gcc
CXXFLAGS :=  -g -O0 -Wall $(Preprocessors)
CFLAGS   :=  -g -O0 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/european_call_binominal_tree.cpp$(ObjectSuffix) $(IntermediateDirectory)/european_put_binomial_tree.cpp$(ObjectSuffix) $(IntermediateDirectory)/binomial_tree.cpp$(ObjectSuffix) $(IntermediateDirectory)/american_put_binominal_tree.cpp$(ObjectSuffix) $(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IntermediateDirectory)/american_call_binominal_tree.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


$(IntermediateDirectory)/.d:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/european_call_binominal_tree.cpp$(ObjectSuffix): european_call_binominal_tree.cpp $(IntermediateDirectory)/european_call_binominal_tree.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/zecophy/Desktop/numercialMethod/9821-Fall-2018-Team-5/HW11/HW11_Tree/european_call_binominal_tree.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/european_call_binominal_tree.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/european_call_binominal_tree.cpp$(DependSuffix): european_call_binominal_tree.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/european_call_binominal_tree.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/european_call_binominal_tree.cpp$(DependSuffix) -MM european_call_binominal_tree.cpp

$(IntermediateDirectory)/european_call_binominal_tree.cpp$(PreprocessSuffix): european_call_binominal_tree.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/european_call_binominal_tree.cpp$(PreprocessSuffix) european_call_binominal_tree.cpp

$(IntermediateDirectory)/european_put_binomial_tree.cpp$(ObjectSuffix): european_put_binomial_tree.cpp $(IntermediateDirectory)/european_put_binomial_tree.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/zecophy/Desktop/numercialMethod/9821-Fall-2018-Team-5/HW11/HW11_Tree/european_put_binomial_tree.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/european_put_binomial_tree.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/european_put_binomial_tree.cpp$(DependSuffix): european_put_binomial_tree.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/european_put_binomial_tree.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/european_put_binomial_tree.cpp$(DependSuffix) -MM european_put_binomial_tree.cpp

$(IntermediateDirectory)/european_put_binomial_tree.cpp$(PreprocessSuffix): european_put_binomial_tree.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/european_put_binomial_tree.cpp$(PreprocessSuffix) european_put_binomial_tree.cpp

$(IntermediateDirectory)/binomial_tree.cpp$(ObjectSuffix): binomial_tree.cpp $(IntermediateDirectory)/binomial_tree.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/zecophy/Desktop/numercialMethod/9821-Fall-2018-Team-5/HW11/HW11_Tree/binomial_tree.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/binomial_tree.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/binomial_tree.cpp$(DependSuffix): binomial_tree.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/binomial_tree.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/binomial_tree.cpp$(DependSuffix) -MM binomial_tree.cpp

$(IntermediateDirectory)/binomial_tree.cpp$(PreprocessSuffix): binomial_tree.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/binomial_tree.cpp$(PreprocessSuffix) binomial_tree.cpp

$(IntermediateDirectory)/american_put_binominal_tree.cpp$(ObjectSuffix): american_put_binominal_tree.cpp $(IntermediateDirectory)/american_put_binominal_tree.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/zecophy/Desktop/numercialMethod/9821-Fall-2018-Team-5/HW11/HW11_Tree/american_put_binominal_tree.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/american_put_binominal_tree.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/american_put_binominal_tree.cpp$(DependSuffix): american_put_binominal_tree.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/american_put_binominal_tree.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/american_put_binominal_tree.cpp$(DependSuffix) -MM american_put_binominal_tree.cpp

$(IntermediateDirectory)/american_put_binominal_tree.cpp$(PreprocessSuffix): american_put_binominal_tree.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/american_put_binominal_tree.cpp$(PreprocessSuffix) american_put_binominal_tree.cpp

$(IntermediateDirectory)/main.cpp$(ObjectSuffix): main.cpp $(IntermediateDirectory)/main.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/zecophy/Desktop/numercialMethod/9821-Fall-2018-Team-5/HW11/HW11_Tree/main.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/main.cpp$(DependSuffix): main.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/main.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/main.cpp$(DependSuffix) -MM main.cpp

$(IntermediateDirectory)/main.cpp$(PreprocessSuffix): main.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/main.cpp$(PreprocessSuffix) main.cpp

$(IntermediateDirectory)/american_call_binominal_tree.cpp$(ObjectSuffix): american_call_binominal_tree.cpp $(IntermediateDirectory)/american_call_binominal_tree.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/zecophy/Desktop/numercialMethod/9821-Fall-2018-Team-5/HW11/HW11_Tree/american_call_binominal_tree.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/american_call_binominal_tree.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/american_call_binominal_tree.cpp$(DependSuffix): american_call_binominal_tree.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/american_call_binominal_tree.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/american_call_binominal_tree.cpp$(DependSuffix) -MM american_call_binominal_tree.cpp

$(IntermediateDirectory)/american_call_binominal_tree.cpp$(PreprocessSuffix): american_call_binominal_tree.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/american_call_binominal_tree.cpp$(PreprocessSuffix) american_call_binominal_tree.cpp


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


