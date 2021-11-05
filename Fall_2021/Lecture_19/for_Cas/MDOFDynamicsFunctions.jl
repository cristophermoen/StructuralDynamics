using LinearAlgebra
using Dierckx
using DiffEqOperators
using DifferentialEquations

function CalculateMemberLength(MemberDefinitions, NodeGeometry)

    L=zeros(length(MemberDefinitions))

    for i=1:length(MemberDefinitions)

        iNode=MemberDefinitions[i][1]
        jNode=MemberDefinitions[i][2]

        iNodePosition=NodeGeometry[iNode,:]
        jNodePosition=NodeGeometry[jNode,:]

        L[i] =norm(iNodePosition-jNodePosition)

    end

    return L

end


function DefineMemberProperty(MemberDefinitions, MemberProperty, PropertyOrder, PropertyType)

    A=zeros(length(MemberDefinitions))

    for i=1:length(MemberDefinitions)

        A[i] = MemberProperty[MemberDefinitions[i][PropertyOrder]][PropertyType]

    end

    return A

end


function CalculateMemberOrientation(MemberDefinitions, NodeGeometry)

    α=zeros(length(MemberDefinitions))

    for i=1:length(MemberDefinitions)

        iNode=MemberDefinitions[i][1]
        jNode=MemberDefinitions[i][2]

        iNodePosition=NodeGeometry[iNode,:]
        jNodePosition=NodeGeometry[jNode,:]

        MemberVector=jNodePosition-iNodePosition
        HorizontalVector=[1; 0]

        α[i] =acos(dot(MemberVector,HorizontalVector)/(norm(MemberVector)*norm(HorizontalVector)))

    end

    return α

end

function DefineLocalKe(A,E,Ix,L)

    k=zeros(6,6)

    k[1,1]=A*E/L
    k[1,4]=-A*E/L
    k[2,2]=12E*Ix/L^3
    k[2,3]=6E*Ix/L^2
    k[2,5]=-12E*Ix/L^3
    k[2,6]=6E*Ix/L^2
    k[3,3]=4E*Ix/L
    k[3,5]=-6E*Ix/L^2
    k[3,6]=2E*Ix/L
    k[4,4]=A*E/L
    k[5,5]=12*E*Ix/L^3
    k[5,6]=-6E*Ix/L^2
    k[6,6]=-4E*Ix/L

    k[2:6,1]=k[1,2:6]
    k[3:6,2]=k[2,3:6]
    k[4:6,3]=k[3,4:6]
    k[5:6,4]=k[4,5:6]
    k[6,5]=k[5,6]

    return k

end


function DefineLocalKe(A,E,Ix,L)

    k=zeros(6,6)

    k[1,1]=A*E/L
    k[1,4]=-A*E/L
    k[2,2]=12E*Ix/L^3
    k[2,3]=6E*Ix/L^2
    k[2,5]=-12E*Ix/L^3
    k[2,6]=6E*Ix/L^2
    k[3,3]=4E*Ix/L
    k[3,5]=-6E*Ix/L^2
    k[3,6]=2E*Ix/L
    k[4,4]=A*E/L
    k[5,5]=12*E*Ix/L^3
    k[5,6]=-6E*Ix/L^2
    k[6,6]=4E*Ix/L

    k[2:6,1]=k[1,2:6]
    k[3:6,2]=k[2,3:6]
    k[4:6,3]=k[3,4:6]
    k[5:6,4]=k[4,5:6]
    k[6,5]=k[5,6]

    return k

end



function DefineRotationMatrix(α)

    θ=zeros(6,6)

    θ[1,1]=cos(α)
    θ[1,2]=sin(α)
    θ[2,1]=-sin(α)
    θ[2,2]=cos(α)
    θ[3,3]=1.0

    θ[4:6,4:6]=θ[1:3,1:3]

    return θ

end

function LocalToGlobal(Local, θ)

    Global=transpose(θ)*Local*θ

end


function DefineLocalMe(A,L,ρ)

    m=zeros(6,6)

    m[1,1]=140
    m[1,4]=70
    m[2,2]=156
    m[2,3]=22L
    m[2,5]=54
    m[2,6]=-13L
    m[3,3]=4L^2
    m[3,5]=13L
    m[3,6]=-3L^2
    m[4,4]=140
    m[5,5]=156
    m[5,6]=-22L
    m[6,6]=4L^2

    m[2:6,1]=m[1,2:6]
    m[3:6,2]=m[2,3:6]
    m[4:6,3]=m[3,4:6]
    m[5:6,4]=m[4,5:6]
    m[6,5]=m[5,6]

    m=(ρ*A*L/420)*m

    return m

end

function DefineFixedDOF(Supports)

    GlobalFixedDOF=zeros(Int64,length(Supports))

    for i=1:length(Supports)

        NodeNumber=Supports[i][1]

        DOFStart=3*(NodeNumber-1)

        GlobalFixedDOF[i]=DOFStart+Supports[i][2]

    end

    return GlobalFixedDOF

end


function AssembleGlobalMatrix(NodeGeometry,MemberDefinitions,eGlobal)

    NumberOfNodes=size(NodeGeometry)[1]

    Global=zeros(NumberOfNodes*3,NumberOfNodes*3)

    for i=1:length(MemberDefinitions)

        iNode=MemberDefinitions[i][1]
        jNode=MemberDefinitions[i][2]

        iNodeDOF=[1;2;3] .+(iNode-1)*3
        jNodeDOF=[1;2;3] .+(jNode-1)*3

        GlobalElementDOF=[iNodeDOF;jNodeDOF]

        Global[GlobalElementDOF,GlobalElementDOF] = Global[GlobalElementDOF,GlobalElementDOF] + eGlobal[i]

    end

    return Global

end

function DefineFreeDOF(Supports)

    FixedDOF = DefineFixedDOF(Supports)

    AllDOF=1:size(NodeGeometry)[1]*3

    FreeDOF=setdiff(AllDOF,FixedDOF)

    return FreeDOF

end



function CalculateGlobalK(MemberDefinitions, SectionProperties, MaterialProperties, NodeGeometry, Supports)

    #calculate member length
    L = CalculateMemberLength(MemberDefinitions, NodeGeometry)

    #calculate member orientation
    α = CalculateMemberOrientation(MemberDefinitions, NodeGeometry)

    #define section properties for stiffness matrix
    A = DefineMemberProperty(MemberDefinitions, SectionProperties, 3, 1)
    E = DefineMemberProperty(MemberDefinitions, MaterialProperties, 4, 1)
    Ix = DefineMemberProperty(MemberDefinitions, SectionProperties, 3, 2)

    #calculate local stiffness matrix for each element
    LocalKe = DefineLocalKe.(A,E,Ix,L)

    #define rotation matrix for each element
    θ = DefineRotationMatrix.(α)

    #rotate local stiffness matrix into global coordinates
    GlobalKe = LocalToGlobal.(LocalKe, θ)

    #assemble global stiffness matrix
    GlobalK = AssembleGlobalMatrix(NodeGeometry,MemberDefinitions,GlobalKe)

    #partition global stiffness matrix, only free DOF
    FreeDOF = DefineFreeDOF(Supports)
    K=GlobalK[FreeDOF,FreeDOF]

    return K

end


function CalculateGlobalM(MemberDefinitions, SectionProperties, MaterialProperties, NodeGeometry, Supports)

    #define section and material properties
    A = DefineMemberProperty(MemberDefinitions, SectionProperties, 3, 1)
    ρ = DefineMemberProperty(MemberDefinitions, MaterialProperties, 4, 2)

    #calculate member length
    L = CalculateMemberLength(MemberDefinitions, NodeGeometry)

    #calculate member orientation
    α = CalculateMemberOrientation(MemberDefinitions, NodeGeometry)

    #define local mass matrix
    LocalMe = DefineLocalMe.(A,L,ρ)

    #define rotation matrix for each element
    θ = DefineRotationMatrix.(α)

    #rotate local mass matrix into global coordinates
    GlobalMe = LocalToGlobal.(LocalMe, θ)

    #assemble global mass matrix
    GlobalM = AssembleGlobalMatrix(NodeGeometry,MemberDefinitions,GlobalMe)

    #partition global mass matrix, only free DOF
    FreeDOF = DefineFreeDOF(Supports)
    M=GlobalM[FreeDOF,FreeDOF]

    return M

end