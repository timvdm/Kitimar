#pragma once

#include <Kitimar/Serialize/Serialize.hpp>
#include <Kitimar/Util/Timer.hpp>

namespace Kitimar {

    inline auto chembl_smi_filename()
    {
        return std::string(KITIMAR_DATA_DIR) + "/chembl_32.smi";
    }

    inline auto chembl_serialized_filename(StructMolecule)
    {
        return std::string(KITIMAR_DATA_DIR) + "/chembl_32.StructMolecule";
    }

    inline auto chembl_serialized_filename(StructMoleculeIncident)
    {
        return std::string(KITIMAR_DATA_DIR) + "/chembl_32.StructMoleculeIncident";
    }

    template<typename Layout>
    void serialize_chembl()
    {
        auto filename = chembl_smi_filename();
        std::ofstream ofs(chembl_serialized_filename(Layout{}), std::ios_base::binary | std::ios_base::out);
        assert(ofs);

        ofs.seekp(LayoutSize::size());

        LayoutSize::Type n = 0;
        std::vector<std::byte> data;
        readMolecules(filename, [&n, &ofs, &data] (auto &mol) {
            ++n;
            if (n % 1000 == 0) std::cout << n << std::endl;
            serialize<Layout>(mol, ofs, data);
            return true;
        });

        //std::cout << "# molecules: " << n << '\n'; // 2327928

        ofs.seekp(0);
        ofs.write(reinterpret_cast<const char*>(&n), LayoutSize::size());
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
