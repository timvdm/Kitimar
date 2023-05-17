#include <Kitimar/RDKit/RDKit.hpp>

int main()
{
    static_assert(Kitimar::AtomList<RDKit::ROMol*>);
    static_assert(Kitimar::BondList<RDKit::ROMol*>);
    static_assert(Kitimar::MoleculeGraph<RDKit::ROMol*>);
    static_assert(Kitimar::IncidentBondList<RDKit::ROMol*>);
    static_assert(Kitimar::AdjacentAtomList<RDKit::ROMol*>);
    static_assert(Kitimar::ElementLayer<RDKit::ROMol*>);
    static_assert(Kitimar::IsotopeLayer<RDKit::ROMol*>);
    static_assert(Kitimar::ChargeLayer<RDKit::ROMol*>);
    static_assert(Kitimar::BondOrderLayer<RDKit::ROMol*>);
    static_assert(Kitimar::ImplicitHydrogensLayer<RDKit::ROMol*>);
    static_assert(Kitimar::AromaticLayer<RDKit::ROMol*>);

}
