#pragma once

#include <Kitimar/CTLayout/Molecule.hpp>
#include <Kitimar/Util/Timer.hpp>

namespace Kitimar {

    inline auto chembl_smi_filename()
    {
        return std::string(KITIMAR_DATA_DIR) + "/chembl_32.smi";
    }

    inline auto chembl_serialized_filename(CTLayout::StructMolecule)
    {
        return std::string(KITIMAR_DATA_DIR) + "/chembl_32.StructMolecule";
    }

    inline auto chembl_serialized_filename(CTLayout::StructMoleculeIncident)
    {
        return std::string(KITIMAR_DATA_DIR) + "/chembl_32.StructMoleculeIncident";
    }

    template<typename Layout>
    void serialize_chembl()
    {
        CTLayout::MoleculeSink<Layout> sink{chembl_serialized_filename(Layout{})};
        readMolecules(chembl_smi_filename(), [&sink] (auto &mol) {
            sink.write(mol);
            if (sink.size() % 1000 == 0) std::cout << sink.size() << std::endl;
            return true;
        });
    }

    /*
    template<typename FromLayout, typename ToLayout>
    void serialize_chembl()
    {
        std::vector<std::byte> idata, odata;
        std::ifstream ifs(chembl_serialized_filename(FromLayout{}), std::ios_base::binary | std::ios_base::in);
        std::ofstream ofs(chembl_serialized_filename(ToLayout{}), std::ios_base::binary | std::ios_base::out);
        assert(ofs);
        while (ifs) {
            auto [ok, mol] = deserialize<FromLayout>(ifs, idata);
            if (!ok)
                break;
            serialize<ToLayout>(mol, ofs, odata);
        }

        ofs.write(reinterpret_cast<const char*>(&LayoutEnd), LayoutSize::size());
    }
*/

    template<typename F>
    auto benchmark(F f)
    {
        Timer timer;
        timer.start();
        if constexpr (std::is_same_v<decltype(f()), void>) {
            f();
            return timer.elapsed();
        } else {
            auto result = f();
            return std::make_pair(result, timer.elapsed());
        }
    }

    template<typename F>
    auto benchmark(const std::string &name, F f)
    {
        Timer timer;
        timer.start();
        f();
        auto elapsed = timer.elapsed();
        std::cout << name << ": " << bechmark(std::forward<F>(f)) << std::endl;
    }



} // namespace Kitimar
