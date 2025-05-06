using ITensors, ITensorMPS
using LinearAlgebra

function ITensors.space(::SiteType"Exciton"; conserve_qns=false, conserve_excitations=conserve_qns, qnname_excitations="excitation")
    if conserve_excitations
        [
            QN(qnname_excitations, 0) => 1,
            QN(qnname_excitations, 1) => 1
        ]
    else
        2
    end
end

ITensors.val(::ValName"G", ::SiteType"Exciton") = 1
ITensors.val(::ValName"E", ::SiteType"Exciton") = 2

ITensors.state(::StateName"G", ::SiteType"Exciton") = [1.0, 0.0]
ITensors.state(::StateName"E", ::SiteType"Exciton") = [0.0, 1.0]

ITensors.op(::OpName"G", s::SiteType"Exciton") = kron(ITensors.state(StateName("G"), s), ITensors.state(StateName("G"), s))
ITensors.op(::OpName"E", s::SiteType"Exciton") = kron(ITensors.state(StateName("E"), s), ITensors.state(StateName("E"), s))
ITensors.op(::OpName"G->E", s::SiteType"Exciton") = kron(ITensors.state(StateName("G"), s), ITensors.state(StateName("E"), s))
ITensors.op(::OpName"E->G", s::SiteType"Exciton") = kron(ITensors.state(StateName("E"), s), ITensors.state(StateName("G"), s))
ITensors.op(::OpName"id", ::SiteType"Exciton") = Matrix(1.0I, 2, 2)