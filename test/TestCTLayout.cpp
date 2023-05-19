#include <Kitimar/CTLayout/CTLayout.hpp>

#include <gtest/gtest.h>

#include <cmath>

using namespace Kitimar;
using namespace Kitimar::CTLayout;

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
