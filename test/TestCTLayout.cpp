#include <Kitimar/CTLayout/Value.hpp>
#include <Kitimar/CTLayout/Struct.hpp>
#include <Kitimar/CTLayout/Vector.hpp>
#include <Kitimar/CTLayout/Sink.hpp>
//#include <Kitimar/CTLayout/CTLayout.hpp>
#include <Kitimar/Util/Util.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <numeric>


#define FMT_HEADER_ONLY
#include <fmt/core.h>

using namespace Kitimar;
using namespace Kitimar::CTLayout;







template<typename T, typename ...Ts, typename Ptr = const std::byte*>
void print(ctll::list<T, Ts...>, Ptr data = nullptr, int depth = 0)
{
    print(T{}, data, depth);
    if constexpr (sizeof...(Ts)) {
        if constexpr (isFixedSize(T{}))
            print(ctll::list<Ts...>{}, data + sizeOf(T{}), depth);
        else
            print(ctll::list<Ts...>{}, data + sizeOf(T{}, data), depth);
    }
}



template<typename T, typename Ptr = const std::byte*>
void print(T, Ptr data = nullptr, int depth = 0)
{
    constexpr auto format = "{:.<50}  {: <8}";
    constexpr auto format2 = "{: <50}  {: <8}";

    if (!depth)
        fmt::println(format2, "TYPE", "SIZE");


    auto indent = std::string(4 * depth, ' ');

    if constexpr (isValue(T{})) {
        auto name = fmt::format("{}{}<{}>  ", indent, Util::typeName(T{}), Util::typeName(typename T::Type{}));
        fmt::println(format, name, sizeOf(T{}));
    }

    if constexpr (isStruct(T{})) {
        auto name = fmt::format("{}{}", indent, "Struct  ");
        if constexpr (isFixedSize(T{}))
            fmt::println(format, name, sizeOf(T{}));
        else
            fmt::println(format, name, sizeOf(T{}, data));
        print(T::members, data, depth + 1);
    }

    if constexpr (isVector(T{})) {
        auto name = fmt::format("{}{}", indent, "Vector  ");
        fmt::println(format, name, sizeOf(T{}, data));
        print(typename T::Type{}, data, depth + 1);
    }

}






/*

CTLayout

- Types
    - Value
        - fixed size
    - Struct
        - fixed size members -> fixed size
        - variable size members -> store size
    - Array
        - store length
        - fixed size value type -> store stride
        - variable size value type -> offset table


- Special Types
    - EndValue (rename to End?)
    - ArraySize<Array>

- isX (X = Value, Struct, ...)
- isFixedSize
- contains
- sizeOf


- Address
    - offset            std::size_t
    - arraySkip         std::array<Size, NumArrays>
    - arrayStride       std::array<Size, ArrayDepth>


- ArraySkip             add size
- ArrayOffset           add offset
- ArrayStride           add i * stride


Types               Fixed

Struct              N
    Value           Y
    Array           N
        Struct      Y
    Array           N
        Array       N
            Value   Y


*/

TEST(TestCTLayout, Value)
{
    struct IntValue : Value<int> {};
    struct BoolValue : Value<bool> {};

    static_assert(isValue(IntValue{}));
    static_assert(!isValue(Struct<IntValue>{}));

    static_assert(isFixedSize(IntValue{}));

    static_assert(contains(IntValue{}, IntValue{}));
    static_assert(!contains(IntValue{}, BoolValue{}));

    static_assert(sizeOf(IntValue{}) == sizeof(int));
    static_assert(sizeOf(BoolValue{}) == sizeof(bool));

    static_assert(offset(IntValue{}, IntValue{}) == 0);
    static_assert(offset(IntValue{}, BoolValue{}) == sizeOf(IntValue{}));
    static_assert(offset(BoolValue{}, IntValue{}) == sizeOf(BoolValue{}));
}

TEST(TestCTLayout, FixedSizeStruct)
{
    struct IntValue : Value<int> {};
    struct ShortValue : Value<short> {};
    struct BoolValue : Value<bool> {};
    struct S : Struct<IntValue, BoolValue> {};

    static_assert(isStruct(S{}));
    static_assert(!isStruct(IntValue{}));

    static_assert(isFixedSize(S{}));

    static_assert(contains(S{}, S{}));
    static_assert(contains(S{}, IntValue{}));
    static_assert(contains(S{}, BoolValue{}));
    static_assert(!contains(S{}, ShortValue{}));

    static_assert(sizeOf(S{}) == sizeOf(IntValue{}) + sizeOf(BoolValue{}));

    static_assert(offset(S{}, S{}) == 0);
    static_assert(offset(S{}, IntValue{}) == 0);
    static_assert(offset(S{}, BoolValue{}) == sizeOf(IntValue{}));
    static_assert(offset(S{}, ShortValue{}) == sizeOf(S{}));
}

TEST(TestCTLayout, FixedSizeVector)
{
    struct IntValue : Value<int> {};
    struct BoolValue : Value<bool> {};
    struct V : Vector<IntValue> {};

    V::Length l = 5;
    std::vector<std::byte> data(sizeof(V::Length) + l * sizeOf(V::Type{}));
    *reinterpret_cast<V::Length*>(data.data()) = l;

    static_assert(isVector(V{}));
    static_assert(!isVector(IntValue{}));

    static_assert(!isFixedSize(V{}));

    static_assert(contains(V{}, V{}));
    static_assert(contains(V{}, IntValue{}));
    static_assert(!contains(V{}, BoolValue{}));

    auto source = PtrSource{data.data()};

    EXPECT_EQ(sizeOf(V{}, source), data.size());

    static_assert(offset(V{}, V{}) == 0);
    EXPECT_EQ(offset(V{}, IntValue{}, source, 0), sizeof(V::Length));
    EXPECT_EQ(offset(V{}, IntValue{}, source, 1), sizeof(V::Length) + sizeOf(IntValue{}));
    EXPECT_EQ(offset(V{}, IntValue{}, source, 2), sizeof(V::Length) + 2 * sizeOf(IntValue{}));
    EXPECT_EQ(offset(V{}, IntValue{}, source, 3), sizeof(V::Length) + 3 * sizeOf(IntValue{}));
    EXPECT_EQ(offset(V{}, IntValue{}, source, 4), sizeof(V::Length) + 4 * sizeOf(IntValue{}));
    EXPECT_EQ(offset(V{}, BoolValue{}, source, 0), sizeOf(V{}, source));
}

TEST(TestCTLayout, VariableSizeStruct)
{
    struct IntValue : Value<int> {};
    struct ShortValue : Value<short> {};
    struct BoolValue : Value<bool> {};
    struct DoubleValue : Value<double> {};
    struct V : Vector<ShortValue> {};
    struct S : Struct<IntValue, V, BoolValue> {};

    // 4 + 4 + 3 * 2 + 1 = 15
    V::Length l = 3;
    std::vector<std::byte> data(sizeOf(IntValue{}) + sizeof(V::Length) + l * sizeOf(V::Type{}) + sizeOf(BoolValue{}));
    *reinterpret_cast<V::Length*>(data.data() + sizeOf(IntValue{})) = l;

    static_assert(!isFixedSize(S{}));

    static_assert(contains(S{}, S{}));
    static_assert(contains(S{}, IntValue{}));
    static_assert(contains(S{}, V{}));
    static_assert(contains(S{}, ShortValue{}));
    static_assert(contains(S{}, BoolValue{}));
    static_assert(!contains(S{}, DoubleValue{}));

    auto source = PtrSource{data.data()};

    std::cout << std::endl;
    print(S{}, source);
    std::cout << std::endl;

    EXPECT_EQ(sizeOf(S{}, source), data.size());

    EXPECT_EQ(offset(S{}, S{}), 0);
    EXPECT_EQ(offset(S{}, IntValue{}, source), 0);
    EXPECT_EQ(offset(S{}, V{}, source), sizeOf(IntValue{}));
    EXPECT_EQ(offset(S{}, ShortValue{}, source, 0), sizeOf(IntValue{}) + sizeof(V::Length));
    EXPECT_EQ(offset(S{}, ShortValue{}, source, 1), sizeOf(IntValue{}) + sizeof(V::Length) + sizeOf(V::Type{}));
    EXPECT_EQ(offset(S{}, ShortValue{}, source, 2), sizeOf(IntValue{}) + sizeof(V::Length) + 2 * sizeOf(V::Type{}));
    EXPECT_EQ(offset(S{}, BoolValue{}, source), sizeOf(IntValue{}) + sizeof(V::Length) + 3 * sizeOf(V::Type{}));
    EXPECT_EQ(offset(S{}, DoubleValue{}, source), sizeOf(S{}, source));

}

TEST(TestCTLayout, VariableSizeVector)
{
    struct ShortValue : Value<short> {};

    struct V2 : Vector<ShortValue> {};
    struct V1 : Vector<V2> {};

    using Length = V1::Length;

    // 4 + 3 * 8 + 3 * 4 + (1 + 2 + 1) * 2 = 48
    Length l1 = 3;
    Length l2[] = {1, 2, 1};
    std::vector<std::byte> data(sizeof(Length) +
                                l1 * sizeof(std::size_t) +
                                l1 * sizeof(Length) +
                                std::accumulate(l2, l2 + l1, 0) * sizeOf(V2::type));
    // V1 length
    *reinterpret_cast<Length*>(data.data()) = l1;
    // V1 offsets
    auto V2_0 = sizeof(Length) + l1 * sizeof(std::size_t);
    auto V2_1 = V2_0 + sizeof(Length) + l2[0] * sizeOf(V2::type);
    auto V2_2 = V2_1 + sizeof(Length) + l2[1] * sizeOf(V2::type);
    *reinterpret_cast<std::size_t*>(data.data() + sizeof(Length)                          ) = 48;
    *reinterpret_cast<std::size_t*>(data.data() + sizeof(Length) +     sizeof(std::size_t)) = V2_1;
    *reinterpret_cast<std::size_t*>(data.data() + sizeof(Length) + 2 * sizeof(std::size_t)) = V2_2;

    *reinterpret_cast<Length*>(data.data() + V2_0) = l2[0];
    *reinterpret_cast<Length*>(data.data() + V2_1) = l2[1];
    *reinterpret_cast<Length*>(data.data() + V2_2) = l2[2];

    static_assert(isVector(V1{}));

    static_assert(!isFixedSize(V1{}));

    static_assert(contains(V1{}, V1{}));
    static_assert(contains(V1{}, V2{}));
    static_assert(contains(V1{}, ShortValue{}));
    static_assert(!contains(V1{}, Value<bool>{}));

    auto source = PtrSource{data.data()};

    EXPECT_EQ(sizeOf(V1{}, source), data.size());

    static_assert(offset(V1{}, V1{}) == 0);
    EXPECT_EQ(offset(V1{}, V2{}, source, 0), V2_0);
    EXPECT_EQ(offset(V1{}, V2{}, source, 1), V2_1);
    EXPECT_EQ(offset(V1{}, V2{}, source, 2), V2_2);
    EXPECT_EQ(offset(V1{}, ShortValue{}, source, 0), V2_0 + sizeof(Length));
    EXPECT_EQ(offset(V1{}, ShortValue{}, source, 1, 0), V2_1 + sizeof(Length));
    EXPECT_EQ(offset(V1{}, ShortValue{}, source, 1, 1), V2_1 + sizeof(Length) + sizeOf(V2::type));
    EXPECT_EQ(offset(V1{}, ShortValue{}, source, 2), V2_2 + sizeof(Length));
    EXPECT_EQ(offset(V1{}, Value<bool>{}, source), sizeOf(V1{}, source));
}


TEST(TestCTLayout, ValueObject)
{
    struct IntValue : Value<int> {};
    std::vector<std::byte> data(sizeOf(IntValue{}));
    auto source = PtrSink{data.data()};
    auto obj = toObject(IntValue{}, source);
    obj.set(42);
    EXPECT_EQ(obj.get(), 42);
}



TEST(TestCTLayout, FixedSizeStructObject)
{
    struct IntValue : Value<int> {};
    struct ShortValue : Value<short> {};
    struct S : Struct<IntValue, ShortValue> {};
    std::vector<std::byte> data(sizeOf(S{}));
    auto source = PtrSink{data.data()};

    auto obj = toObject(S{}, source);

    auto intObj = obj.get(IntValue{});
    auto shortObj = obj.get(ShortValue{});
    //EXPECT_EQ(intObj.data(), obj.data());
    //EXPECT_EQ(shortObj.data(), obj.data() + sizeOf(IntValue{}));

    intObj.set(42);
    shortObj.set(24);
    EXPECT_EQ(intObj.get(), 42);
    EXPECT_EQ(shortObj.get(), 24);
}


TEST(TestCTLayout, FixedSizeVectorObject)
{
    struct IntValue : Value<int> {};
    struct V : Vector<IntValue> {};

    V::Length l = 3;
    std::vector<std::byte> data(sizeof(V::Length) + l * sizeOf(V::Type{}));
    *reinterpret_cast<V::Length*>(data.data()) = l;
    *reinterpret_cast<IntValue::Type*>(data.data() + sizeof(V::Length)) = 1;
    *reinterpret_cast<IntValue::Type*>(data.data() + sizeof(V::Length) + sizeof(IntValue::Type)) = 2;
    *reinterpret_cast<IntValue::Type*>(data.data() + sizeof(V::Length) + 2 * sizeof(IntValue::Type)) = 3;

    auto source = PtrSink{data.data()};

    auto v = toObject(V{}, source);
    EXPECT_EQ(v.length(), l);
    EXPECT_EQ(v.at(0).get(), 1);
    EXPECT_EQ(v.at(1).get(), 2);
    EXPECT_EQ(v.at(2).get(), 3);

    v.at(0).set(42);
    v.at(1).set(24);
    v.at(2).set(48);

    EXPECT_EQ(v.at(0).get(), 42);
    EXPECT_EQ(v.at(1).get(), 24);
    EXPECT_EQ(v.at(2).get(), 48);
}


TEST(TestCTLayout, VariableSizeStructObject)
{
    struct IntValue : Value<int> {};
    struct ShortValue : Value<short> {};
    struct BoolValue : Value<bool> {};
    struct V : Vector<ShortValue> {};
    struct S : Struct<IntValue, V, BoolValue> {};

    // 4 + 4 + 3 * 2 + 1 = 15
    V::Length l = 3;
    std::vector<std::byte> data(sizeOf(IntValue{}) + sizeof(V::Length) + l * sizeOf(V::Type{}) + sizeOf(BoolValue{}));
    auto ptr = data.data();
    *reinterpret_cast<IntValue::Type*>(ptr) = 1;
    ptr += sizeOf(IntValue{});
    *reinterpret_cast<V::Length*>(ptr) = l;
    ptr += sizeof(V::Length);
    *reinterpret_cast<ShortValue::Type*>(ptr) = 2;
    ptr += sizeOf(ShortValue{});
    *reinterpret_cast<ShortValue::Type*>(ptr) = 3;
    ptr += sizeOf(ShortValue{});
    *reinterpret_cast<ShortValue::Type*>(ptr) = 4;
    ptr += sizeOf(ShortValue{});
    *reinterpret_cast<BoolValue::Type*>(ptr) = true;

    auto source = PtrSink{data.data()};

    auto obj = toObject(S{}, source);
    auto v = obj.get(V{});
    EXPECT_EQ(obj.get(IntValue{}).get(), 1);
    EXPECT_EQ(v.length(), l);
    EXPECT_EQ(v.at(0).get(), 2);
    EXPECT_EQ(v.at(1).get(), 3);
    EXPECT_EQ(v.at(2).get(), 4);
    EXPECT_EQ(obj.get(BoolValue{}).get(), true);


    obj.get(IntValue{}).set(1000);
    v.at(0).set(2000);
    v.at(1).set(3000);
    v.at(2).set(4000);
    obj.get(BoolValue{}).set(false);

    EXPECT_EQ(obj.get(IntValue{}).get(), 1000);
    EXPECT_EQ(v.length(), l);
    EXPECT_EQ(v.at(0).get(), 2000);
    EXPECT_EQ(v.at(1).get(), 3000);
    EXPECT_EQ(v.at(2).get(), 4000);
    EXPECT_EQ(obj.get(BoolValue{}).get(), false);
}

TEST(TestCTLayout, VariableSizeVectorObject)
{
    struct ShortValue : Value<short> {};

    struct V2 : Vector<ShortValue> {};
    struct V1 : Vector<V2> {};

    using Length = V1::Length;

    // 4 + 3 * 8 + 3 * 4 + (1 + 2 + 1) * 4 = 56
    Length l1 = 3;
    Length l2[] = {1, 2, 1};
    std::vector<std::byte> data(sizeof(Length) +
                                l1 * sizeof(std::size_t) +
                                l1 * sizeof(Length) +
                                std::accumulate(l2, l2 + l1, 0) * sizeOf(V2::type));
    // V1 length
    *reinterpret_cast<Length*>(data.data()) = l1;
    // V1 offsets
    auto V2_0 = sizeof(Length) + l1 * sizeof(std::size_t);
    auto V2_1 = V2_0 + sizeof(Length) + l2[0] * sizeOf(V2::type);
    auto V2_2 = V2_1 + sizeof(Length) + l2[1] * sizeOf(V2::type);
    *reinterpret_cast<std::size_t*>(data.data() + sizeof(Length)                          ) = V2_0; // FIXME: store sizeOF here...
    *reinterpret_cast<std::size_t*>(data.data() + sizeof(Length) +     sizeof(std::size_t)) = V2_1;
    *reinterpret_cast<std::size_t*>(data.data() + sizeof(Length) + 2 * sizeof(std::size_t)) = V2_2;

    *reinterpret_cast<Length*>(data.data() + V2_0) = l2[0];
    *reinterpret_cast<Length*>(data.data() + V2_1) = l2[1];
    *reinterpret_cast<Length*>(data.data() + V2_2) = l2[2];

    *reinterpret_cast<ShortValue::Type*>(data.data() + V2_0 + sizeof(Length)) = 10;
    *reinterpret_cast<ShortValue::Type*>(data.data() + V2_1 + sizeof(Length)) = 20;
    *reinterpret_cast<ShortValue::Type*>(data.data() + V2_1 + sizeof(Length) + sizeOf(V2::Type{})) = 30;
    *reinterpret_cast<ShortValue::Type*>(data.data() + V2_2 + sizeof(Length)) = 40;

    auto source = PtrSink{data.data()};

    auto v1 = toObject(V1{}, source);

    EXPECT_EQ(v1.length(), l1);
    EXPECT_EQ(v1.at(0).length(), l2[0]);
    EXPECT_EQ(v1.at(1).length(), l2[1]);
    EXPECT_EQ(v1.at(2).length(), l2[2]);

    EXPECT_EQ(v1.at(0).at(0).get(), 10);
    EXPECT_EQ(v1.at(1).at(0).get(), 20);
    EXPECT_EQ(v1.at(1).at(1).get(), 30);
    EXPECT_EQ(v1.at(2).at(0).get(), 40);


    v1.at(0).at(0).set(1000);
    v1.at(1).at(0).set(2000);
    v1.at(1).at(1).set(3000);
    v1.at(2).at(0).set(4000);

    EXPECT_EQ(v1.at(0).at(0).get(), 1000);
    EXPECT_EQ(v1.at(1).at(0).get(), 2000);
    EXPECT_EQ(v1.at(1).at(1).get(), 3000);
    EXPECT_EQ(v1.at(2).at(0).get(), 4000);
}


TEST(TestCTLayout, ValueWriter)
{
    struct IntValue : Value<int> {};

    std::vector<std::byte> data;
    StlVectorSink sink{data};

    ValueWriter writer{IntValue{}, sink};
    writer.write(42);
    ASSERT_EQ(data.size(), sizeOf(IntValue{}));

    auto source = PtrSource{data.data()};
    auto obj = toObject(IntValue{}, source);
    ASSERT_EQ(obj.get(), 42);
}

TEST(TestCTLayout, FileStreamSink)
{
    auto path = "test_FileStreamSink.bin";
    FileStreamSink sink{path};
    sink.write(2, sizeof(int));
    sink.write(1);
    sink.close();

    auto data = Util::readFileData(path);
    EXPECT_EQ(data.size(), 2 * sizeof(int));
    EXPECT_EQ(*reinterpret_cast<int*>(data.data()), 1);
    EXPECT_EQ(*reinterpret_cast<int*>(data.data() + sizeof(int)), 2);
}


TEST(TestCTLayout, StructWriter)
{
    struct A : Value<int> {};
    struct B : Value<short> {};
    struct S : Struct<A, B> {};

    std::vector<std::byte> data;
    StlVectorSink sink{data};

    StructWriter writer{S{}, sink};
    writer.get(A{}).write(1);
    writer.get(B{}).write(2);    
    EXPECT_EQ(data.size(), sizeOf(S{}));

    auto source = PtrSource{data.data()};
    auto obj = toObject(S{}, source);
    EXPECT_EQ(obj.get(A{}).get(), 1);
    EXPECT_EQ(obj.get(B{}).get(), 2);
}

TEST(TestCTLayout, FixedSizeVectorWriter)
{
    struct T : Value<int> {};
    struct V : Vector<T> {};

    std::vector<std::byte> data;
    StlVectorSink sink{data};

    VectorWriter writer{V{}, sink};
    writer.setLength(3);
    EXPECT_EQ(data.size(), sizeof(V::Length));
    EXPECT_EQ(writer.length(), 3);

    writer.at(0).write(10);
    writer.at(1).write(20);
    writer.at(2).write(30);

    auto source = PtrSource{data.data()};
    auto v = toObject(V{}, source);
    EXPECT_EQ(v.length(), 3);
    EXPECT_EQ(v.at(0).get(), 10);
    EXPECT_EQ(v.at(1).get(), 20);
    EXPECT_EQ(v.at(2).get(), 30);
}


/*

                offset  size    value
  length        0       4       3       OK
  sizeOf        4       8       56
  offset1       12      8       36      OK
  offset2       20      8       48
  length0       28      4       1       OK
  value00       32      4       10      OK
  length1       36      4       2       OK
  value10       40      4       20      OK
  value11       44      4       30      OK
  length2       48      4       1       OK
  value20       52      4       40      OK


*/

TEST(TestCTLayout, VariableSizeVectorWriter)
{
    struct T : Value<int> {};
    struct W : Vector<T> {};
    struct V : Vector<W> {};

    std::vector<std::byte> data;
    StlVectorSink sink{data};

    VectorWriter writer{V{}, sink};
    writer.setLength(3);
    EXPECT_EQ(data.size(), sizeof(V::Length));
    EXPECT_EQ(writer.length(), 3);

    auto v0 = writer.at(0);
    v0.setLength(1);
    EXPECT_EQ(v0.length(), 1);
    v0.at(0).write(10);



    auto v1 = writer.at(1);
    v1.setLength(2);
    EXPECT_EQ(v1.length(), 2);    
    v1.at(0).write(20);
    v1.at(1).write(30);    

    auto v2 = writer.at(2);
    v2.setLength(1);
    EXPECT_EQ(v2.length(), 1);
    v2.at(0).write(40);

    EXPECT_EQ(writer.length(), 3);

    auto source = PtrSource{data.data()};
    auto v = toObject(V{}, source);
    EXPECT_EQ(v.length(), 3);
    EXPECT_EQ(v.at(0).length(), 1);
    EXPECT_EQ(v.at(1).length(), 2);
    EXPECT_EQ(v.at(2).length(), 1);
    EXPECT_EQ(v.at(0).at(0).get(), 10);
    EXPECT_EQ(v.at(1).at(0).get(), 20);
    EXPECT_EQ(v.at(1).at(1).get(), 30);
    EXPECT_EQ(v.at(2).at(0).get(), 40);    
    EXPECT_EQ(sizeOf(V{}, sink), 56);
}




















/*



// Value
struct BoolValue : Value<bool> {};
struct CharValue : Value<char> {};
struct ShortValue : Value<short> {};
struct IntValue : Value<int> {};
struct FloatValue : Value<float> {};
struct DoubleValue : Value<double> {};

// Struct
struct EmptyStruct : Struct<> {};
struct ValueStruct : Struct<CharValue, ShortValue, IntValue> {};
struct RecursiveStruct : Struct<FloatValue, ValueStruct> {};

// Array
struct ArrayValue : Value<int> {};
struct NumValues : ArraySize {};
struct ValueArray : Array<ArrayValue, NumValues> {};
struct NumStructs : ArraySize {};
struct StructArray : Array<ValueStruct, NumStructs> {};


TEST(TestCTLayout, IsType)
{
    // isValue(T)
    static_assert( isValue(IntValue()));
    static_assert(!isValue(ValueStruct()));
    static_assert(!isValue(ValueArray()));

    // isStruct(T)
    static_assert(!isStruct(IntValue()));
    static_assert( isStruct(ValueStruct()));
    static_assert(!isStruct(ValueArray()));

    // isArray(T)
    static_assert(!isArray(IntValue()));
    static_assert(!isArray(ValueStruct()));
    static_assert( isArray(ValueArray()));
}

TEST(TestCTLayout, Size)
{
    // Value<T>::size()
    static_assert(CharValue::size() == sizeof(char));
    static_assert(ShortValue::size() == sizeof(short));
    static_assert(IntValue::size() == sizeof(int));
    // Struct<Value<T>, ...>::size()
    static_assert(ValueStruct::size() == sizeof(char) + sizeof(short) + sizeof(int));
    // Struct<Value<T>, Struct<Ts...>, ...>::size()
    static_assert(RecursiveStruct::size() == sizeof(float) + sizeof(char) + sizeof(short) + sizeof(int));
    // Array<Value<T>>::size()
    static_assert(ValueArray::stride() == ArrayValue::size());
    // Array<Struct<T>>::size()
    static_assert(StructArray::stride() == ValueStruct::size());
}

TEST(TestCTLayout, Find)
{
    // find(Value<T>, Value<T>)
    static_assert(detail::find(IntValue(), IntValue()) == 0);
    static_assert(detail::find(ShortValue(), IntValue()) == IntValue::size());
    static_assert(detail::find(IntValue(), ShortValue()) == ShortValue::size());

    // find(Struct<Ts...>, Struct<Ts...>)
    static_assert(detail::find(ValueStruct(), ValueStruct()) == 0);
    static_assert(detail::find(ValueStruct(), RecursiveStruct()) == FloatValue::size());
    static_assert(detail::find(RecursiveStruct(), ValueStruct()) == ValueStruct::size());
    // find(Value<T>, Struct<Value<T>, ...>)
    static_assert(detail::find(CharValue(), ValueStruct()) == 0);
    static_assert(detail::find(ShortValue(), ValueStruct()) == CharValue::size());
    static_assert(detail::find(IntValue(), ValueStruct()) == CharValue::size() + ShortValue::size());
    static_assert(detail::find(FloatValue(), ValueStruct()) == ValueStruct::size());
    // find(Value<T>, Struct<Value<T>, Struct<Ts...>, ...>)
    static_assert(detail::find(FloatValue(), RecursiveStruct()) == 0);
    static_assert(detail::find(CharValue(), RecursiveStruct()) == FloatValue::size());
    static_assert(detail::find(ShortValue(), RecursiveStruct()) == FloatValue::size() + CharValue::size());
    static_assert(detail::find(IntValue(), RecursiveStruct()) == FloatValue::size() + CharValue::size() + ShortValue::size());
    static_assert(detail::find(BoolValue(), RecursiveStruct()) == RecursiveStruct::size());

    // find(Array<T>, Array<T>)
    static_assert(detail::find(ValueArray(), ValueArray()) == 0);
    static_assert(detail::find(ValueArray(), StructArray()) == StructArray::stride());
    static_assert(detail::find(StructArray(), ValueArray()) == ValueArray::stride());

    // find(Value<T>, Array<Value<T>>)
    static_assert(detail::find(ArrayValue(), ValueArray()) == 0);
    static_assert(detail::find(ShortValue(), ValueArray()) == ValueArray::stride());
    // find(Value<T>, Array<Struct<T>>)
    static_assert(detail::find(CharValue(), StructArray()) == 0);
    static_assert(detail::find(ShortValue(), StructArray()) == CharValue::size());
    static_assert(detail::find(IntValue(), StructArray()) == CharValue::size() + ShortValue::size());
    static_assert(detail::find(FloatValue(), StructArray()) == StructArray::stride());

    // findArrays(...)
    static_assert(std::is_same_v<decltype(detail::findArrays(ctll::empty_list())), ctll::empty_list>);
    static_assert(std::is_same_v<decltype(detail::findArrays(ctll::list<ValueArray>())), ctll::list<ValueArray>>);
    static_assert(std::is_same_v<decltype(detail::findArrays(ctll::list<IntValue, ValueArray>())), ctll::list<ValueArray>>);
    static_assert(std::is_same_v<decltype(detail::findArrays(ctll::list<IntValue, ValueArray, ValueStruct, StructArray>())), ctll::list<ValueArray, StructArray>>);
}

auto set_data = [] <typename T> (std::byte *ptr, auto offset, T value) {
    auto p = reinterpret_cast<T*>(ptr + offset);
    *p = value;
};

TEST(TestCTLayout, Layout1)
{
    struct Layout1 : Layout<FloatValue, ValueArray, DoubleValue, StructArray, BoolValue> {};

    std::cout << Layout1::find(FloatValue{}) << std::endl;
    std::cout << Layout1::find(ValueArray{}) << std::endl;
    std::cout << Layout1::find(ArrayValue{}) << std::endl;
    std::cout << Layout1::find(DoubleValue{}) << std::endl;
    std::cout << Layout1::find(StructArray{}) << std::endl;
    std::cout << Layout1::find(ValueStruct{}) << std::endl;
    std::cout << Layout1::find(CharValue{}) << std::endl;
    std::cout << Layout1::find(ShortValue{}) << std::endl;
    std::cout << Layout1::find(IntValue{}) << std::endl;
    std::cout << Layout1::find(BoolValue{}) << std::endl;
    std::cout << Layout1::find(detail::EndValue{}) << std::endl;

    // Layout::header
    static_assert(std::is_same_v<Layout1::Header, ctll::list<LayoutSize, NumValues, NumStructs>>);
    static_assert(Layout1::headerSize == sizeof(std::size_t) + 2 * sizeof(SizeT));

    // Layout::size(sizes)
    static_assert(Layout1::size({0, 0}) == Layout1::headerSize + FloatValue::size() + DoubleValue::size() + BoolValue::size());
    static_assert(Layout1::size({1, 0}) == Layout1::headerSize + FloatValue::size() + DoubleValue::size() + BoolValue::size() + ValueArray::stride());
    static_assert(Layout1::size({0, 1}) == Layout1::headerSize + FloatValue::size() + DoubleValue::size() + BoolValue::size() + StructArray::stride());
    static_assert(Layout1::size({1, 1}) == Layout1::headerSize + FloatValue::size() + DoubleValue::size() + BoolValue::size() + ValueArray::stride() + StructArray::stride());
    static_assert(Layout1::size({2, 3}) == Layout1::headerSize + FloatValue::size() + DoubleValue::size() + BoolValue::size() + 2 * ValueArray::stride() + 3 * StructArray::stride());

    constexpr SizeT sizes[] = {2, 2};

    // Layout::offset(sizes)
    static_assert(Layout1::offset<FloatValue>(sizes) == Layout1::headerSize);
    static_assert(Layout1::offset<ValueArray>(sizes) == Layout1::headerSize + FloatValue::size());
    static_assert(Layout1::offset<DoubleValue>(sizes) == Layout1::headerSize + FloatValue::size() + 2 * ValueArray::stride() );
    static_assert(Layout1::offset<StructArray>(sizes) == Layout1::headerSize + FloatValue::size() + 2 * ValueArray::stride() + DoubleValue::size());
    static_assert(Layout1::offset<BoolValue>(sizes) == Layout1::headerSize + FloatValue::size() + 2 * ValueArray::stride() + DoubleValue::size() + 2 * StructArray::stride());

    // Layout::offset(sizes, index)
    static_assert(Layout1::offset<ArrayValue>(sizes, 0) == Layout1::headerSize + FloatValue::size());
    static_assert(Layout1::offset<ArrayValue>(sizes, 1) == Layout1::headerSize + FloatValue::size() + ArrayValue::size());
    static_assert(Layout1::offset<CharValue>(sizes, 0) == Layout1::headerSize + FloatValue::size() + 2 * ArrayValue::size() + DoubleValue::size());
    static_assert(Layout1::offset<CharValue>(sizes, 1) == Layout1::headerSize + FloatValue::size() + 2 * ArrayValue::size() + DoubleValue::size() + ValueStruct::size());
    static_assert(Layout1::offset<ShortValue>(sizes, 0) == Layout1::headerSize + FloatValue::size() + 2 * ArrayValue::size() + DoubleValue::size() + CharValue::size());
    static_assert(Layout1::offset<ShortValue>(sizes, 1) == Layout1::headerSize + FloatValue::size() + 2 * ArrayValue::size() + DoubleValue::size() + ValueStruct::size() + CharValue::size());
    static_assert(Layout1::offset<IntValue>(sizes, 0) == Layout1::headerSize + FloatValue::size() + 2 * ArrayValue::size() + DoubleValue::size() + CharValue::size() + ShortValue::size());
    static_assert(Layout1::offset<IntValue>(sizes, 1) == Layout1::headerSize + FloatValue::size() + 2 * ArrayValue::size() + DoubleValue::size() + ValueStruct::size() + CharValue::size() + ShortValue::size());


    //
    // Object<Layout1>
    //

    std::array<std::byte, Layout1::size(sizes)> data1 = {};
    auto ptr1 = data1.data();

    auto is_approx = [] (auto a, auto b) { return std::fabs(a - b) < 10e-7; };

    set_data(ptr1, LayoutSize::size(), sizes[0]);
    set_data(ptr1, LayoutSize::size() + sizeof(SizeT), sizes[1]);
    set_data(ptr1, Layout1::offset<FloatValue>(sizes), float{1.2});
    set_data(ptr1, Layout1::offset<ArrayValue>(sizes, 0), int{3});
    set_data(ptr1, Layout1::offset<ArrayValue>(sizes, 1), int{4});
    set_data(ptr1, Layout1::offset<DoubleValue>(sizes), double{3.4});
    set_data(ptr1, Layout1::offset<CharValue>(sizes, 0), char{'a'});
    set_data(ptr1, Layout1::offset<CharValue>(sizes, 1), char{'b'});
    set_data(ptr1, Layout1::offset<ShortValue>(sizes, 0), short{5});
    set_data(ptr1, Layout1::offset<ShortValue>(sizes, 1), short{6});
    set_data(ptr1, Layout1::offset<IntValue>(sizes, 0), int{7});
    set_data(ptr1, Layout1::offset<IntValue>(sizes, 1), int{8});
    set_data(ptr1, Layout1::offset<BoolValue>(sizes), bool{true});

    auto obj1 = Object<Layout1>{ptr1};

    assert(obj1.size() == Layout1::size(sizes));

    assert(is_approx(obj1.value<FloatValue>(), 1.2));
    assert(obj1.value<ArrayValue>(0) == 3);
    assert(obj1.value<ArrayValue>(1) == 4);
    assert(is_approx(obj1.value<DoubleValue>(), 3.4));
    assert(obj1.value<CharValue>(0) == 'a');
    assert(obj1.value<CharValue>(1) == 'b');
    assert(obj1.value<ShortValue>(0) == 5);
    assert(obj1.value<ShortValue>(1) == 6);
    assert(obj1.value<IntValue>(0) == 7);
    assert(obj1.value<IntValue>(1) == 8);
    assert(obj1.value<BoolValue>() == true);
}


TEST(TestCTLayout, Layout2)
{
    //
    // Array::n > 1
    //
    struct ArrayN1 : Array<IntValue, NumValues> {};
    struct ArrayN2 : Array<ShortValue, NumValues, 2> {};
    struct Layout2 : Layout<ArrayN1, ArrayN2> {};

    static_assert(ArrayN1::stride() == IntValue::size());
    static_assert(ArrayN2::stride() == ShortValue::size());

    static_assert(std::is_same_v<Layout2::Header, ctll::list<LayoutSize, NumValues>>);
    static_assert(Layout2::headerSize == sizeof(std::size_t) + sizeof(SizeT));

    constexpr SizeT sizes[] = {2};

    static_assert(Layout2::size(sizes) == Layout2::headerSize + 2 * IntValue::size() + 4 * ShortValue::size());

    static_assert(Layout2::offset<IntValue>(sizes, 0) == Layout2::headerSize);
    static_assert(Layout2::offset<IntValue>(sizes, 1) == Layout2::headerSize + IntValue::size());
    static_assert(Layout2::offset<ShortValue>(sizes, 0) == Layout2::headerSize + 2 * IntValue::size());
    static_assert(Layout2::offset<ShortValue>(sizes, 1) == Layout2::headerSize + 2 * IntValue::size() + ShortValue::size());
    static_assert(Layout2::offset<ShortValue>(sizes, 2) == Layout2::headerSize + 2 * IntValue::size() + 2 * ShortValue::size());
    static_assert(Layout2::offset<ShortValue>(sizes, 3) == Layout2::headerSize + 2 * IntValue::size() + 3 * ShortValue::size());

    std::array<std::byte, Layout2::size(sizes)> data2 = {};
    auto ptr2 = data2.data();

    set_data(ptr2, 0, Layout2::size(sizes));
    set_data(ptr2, LayoutSize::size(), std::size_t{sizes[0]});
    set_data(ptr2, Layout2::offset<IntValue>(sizes, 0), int{1});
    set_data(ptr2, Layout2::offset<IntValue>(sizes, 1), int{2});
    set_data(ptr2, Layout2::offset<ShortValue>(sizes, 0), short{3});
    set_data(ptr2, Layout2::offset<ShortValue>(sizes, 1), short{4});
    set_data(ptr2, Layout2::offset<ShortValue>(sizes, 2), short{5});
    set_data(ptr2, Layout2::offset<ShortValue>(sizes, 3), short{6});

    auto obj2 = Object<Layout2>{ptr2};

    assert(obj2.size() == LayoutSize::size() + sizeof(SizeT) + 2 * sizeof(int) + 4 * sizeof(short));
    assert(obj2.size<ArrayN1>() == 2);
    assert(obj2.size<ArrayN2>() == 4);

    assert(obj2.value<IntValue>(0) == 1);
    assert(obj2.value<IntValue>(1) == 2);
    assert(obj2.value<ShortValue>(0) == 3);
    assert(obj2.value<ShortValue>(1) == 4);
    assert(obj2.value<ShortValue>(2) == 5);
    assert(obj2.value<ShortValue>(3) == 6);
}

TEST(TestCTLayout, BitArray)
{
    BitArray<NumValues> bits;

    static_assert(bits.size(0) == 0);
    static_assert(bits.size(5) == 1);
    static_assert(bits.size(8) == 1);
    static_assert(bits.size(9) == 2);
    static_assert(bits.size(16) == 2);
    static_assert(bits.size(17) == 3);
    static_assert(bits.size(20) == 3);




    std::byte data3[2] = {};

    assert(data3[0] == std::byte{0});
    assert(data3[1] == std::byte{0});

    bits.set(data3, 0);
    assert(data3[0] == std::byte{1});
    bits.set(data3, 3);
    assert(data3[0] == std::byte{1+8});
    bits.set(data3, 7);
    assert(data3[0] == std::byte{1+8+128});
    assert(data3[1] == std::byte{0});

    bits.set(data3, 8);
    assert(data3[1] == std::byte{1});
    bits.set(data3, 15);
    assert(data3[1] == std::byte{1+128});
    assert(data3[0] == std::byte{1+8+128});

    bits.unset(data3, 8);
    assert(data3[1] == std::byte{128});
    bits.unset(data3, 3);
    assert(data3[0] == std::byte{1+128});
    assert(data3[1] == std::byte{128});

    assert(bits.get(data3, 0));
    assert(!bits.get(data3, 1));
    assert(!bits.get(data3, 3));
    assert(bits.get(data3, 7));
    assert(!bits.get(data3, 8));
    assert(bits.get(data3, 15));


    struct NumBits: ArraySize {};
    struct BitsLayout : Layout<BitArray<NumBits>, IntValue> {};

    constexpr SizeT numBits = 15;
    static_assert(BitsLayout::size(&numBits) == BitsLayout::headerSize + 2 + IntValue::size());
    static_assert(BitsLayout::offset<BitArray<NumBits>>(&numBits) == BitsLayout::headerSize);
    static_assert(BitsLayout::offset<IntValue>(&numBits) == BitsLayout::headerSize + 2);


//    auto bitObj = Object<BitsLayout>{data3};

//    auto bitArray = BitArray<NumBits>{};
//    assert(bitObj.bit(bitArray, 0));
//    assert(!bitObj.bit(bitArray, 1));
//    assert(!bitObj.bit(bitArray, 3));
//    assert(bitObj.bit(bitArray, 7));
//    assert(!bitObj.bit(bitArray, 8));
//    assert(bitObj.bit(bitArray, 15));
}
*/
