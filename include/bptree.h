#pragma once

#include <algorithm>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

template <class Key, class Value, std::size_t BlockSize = 4096, class Less = std::less<Key>>
class BPTree
{
    struct Node;
    struct Leaf;
    struct Inner;

    template <bool is_const>
    class Iterator;

public:
    using key_type = Key;
    using mapped_type = Value;
    using value_type = std::pair<Key, Value>; // NB: a digression from std::map
    using reference = value_type &;
    using const_reference = const value_type &;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using size_type = std::size_t;

    using iterator = Iterator<false>;
    using const_iterator = Iterator<true>;

    BPTree();
    BPTree(std::initializer_list<std::pair<Key, Value>>);
    BPTree(const BPTree &);
    BPTree(BPTree &&);

    BPTree & operator=(const BPTree &);
    BPTree & operator=(BPTree &&);

    iterator begin();
    const_iterator cbegin() const;
    const_iterator begin() const;
    iterator end();
    const_iterator cend() const;
    const_iterator end() const;

    bool empty() const;
    size_type size() const;
    void clear();

    size_type count(const Key &) const;
    bool contains(const Key &) const;
    std::pair<iterator, iterator> equal_range(const Key &);
    std::pair<const_iterator, const_iterator> equal_range(const Key &) const;
    iterator lower_bound(const Key &);
    const_iterator lower_bound(const Key &) const;
    iterator upper_bound(const Key &);
    const_iterator upper_bound(const Key &) const;
    iterator find(const Key & key);
    const_iterator find(const Key & key) const;

    // 'at' method throws std::out_of_range if there is no such key
    Value & at(const Key &);
    const Value & at(const Key &) const;

    // '[]' operator inserts a new element if there is no such key
    Value & operator[](const Key &);

    std::pair<iterator, bool> insert(const Key &, const Value &); // NB: a digression from std::map
    std::pair<iterator, bool> insert(const Key &, Value &&);      // NB: a digression from std::map
    template <class ForwardIt>
    void insert(ForwardIt, ForwardIt);
    void insert(std::initializer_list<value_type>);
    iterator erase(const_iterator);
    iterator erase(const_iterator, const_iterator);
    size_type erase(const Key &);

    ~BPTree();

private:
    static constexpr size_type m_order = std::max<std::size_t>(2, (BlockSize - sizeof(std::size_t) - sizeof(void *)) / (sizeof(Key) + sizeof(void *)));
    Node * m_root = nullptr;
    size_type m_size = 0;

    static bool less(const Key &, const Key &);

    static Leaf * find_min_leaf(Node *);
    static Leaf * find_max_leaf(Node *);

    static Node * copy(Node *);
    void redirect();

    static void clear(Node *);

    std::pair<Leaf *, size_type> search(const Key &) const;

    template <class K, class V>
    iterator insert(Leaf *, size_type, K &&, V &&);

    void split(Leaf *);
    void split(Inner *);

    static Leaf * find_left(Leaf *);
    static Inner * find_left(Inner *);
    static Leaf * find_right(Leaf *);
    static Inner * find_right(Inner *);

    bool erase(Leaf *, size_type, const Key &);
    static void erase(Inner *, size_type, const Key &);

    static void merge(Leaf *, Leaf *, const Key &);
    static void merge(Inner *, Inner *, const Key &);

    static void update(Node *, const Key &, const Key &);
};

template <class Key, class Value, std::size_t BlockSize, class Less>
struct BPTree<Key, Value, BlockSize, Less>::Node
{
    Inner * parent = nullptr;

    virtual size_type size() const = 0;
    virtual const Key & get_key(size_type) const = 0;
    virtual Node * get_child(size_type) const = 0;

    virtual ~Node() = default;
};

template <class Key, class Value, std::size_t BlockSize, class Less>
struct BPTree<Key, Value, BlockSize, Less>::Leaf : Node
{
    Leaf * right = nullptr;
    std::vector<value_type> data;

    Leaf()
        : data(0)
    {
        data.reserve(m_order);
    }

    Leaf(const Leaf & node)
        : data(node.data)
    {
    }

    size_type size() const override
    {
        return data.size();
    }

    const Key & get_key(size_type ind) const override
    {
        if (ind < size()) {
            return data[ind].first;
        }
        throw std::out_of_range("");
    }

    Node * get_child(size_type) const override
    {
        return nullptr;
    }
};

template <class Key, class Value, std::size_t BlockSize, class Less>
struct BPTree<Key, Value, BlockSize, Less>::Inner : Node
{
    std::vector<key_type> keys;
    std::vector<Node *> children;

    Inner()
        : keys(0)
        , children(0)
    {
        keys.reserve(m_order);
        children.reserve(m_order + 1);
    }

    Inner(const Inner & node)
        : keys(node.keys)
        , children(0)
    {
    }

    size_type size() const override
    {
        return keys.size();
    }

    const Key & get_key(size_type ind) const override
    {
        if (ind < size()) {
            return keys[ind];
        }
        throw std::out_of_range("");
    }

    Node * get_child(size_type ind) const override
    {
        return (ind < children.size()) ? children[ind] : nullptr;
    }
};

template <class Key, class Value, std::size_t BlockSize, class Less>
template <bool is_const>
class BPTree<Key, Value, BlockSize, Less>::Iterator
{
public:
    using difference_type = std::ptrdiff_t;
    using value_type = std::conditional_t<is_const, const std::pair<Key, Value>, std::pair<Key, Value>>;
    using pointer = value_type *;
    using reference = value_type &;
    using iterator_category = std::forward_iterator_tag;

    Iterator() = default;

    template <bool OtherType>
    Iterator(const Iterator<OtherType> & it)
        : m_current(it.m_current)
        , m_pos(it.m_pos)
    {
    }

    template <bool OtherType>
    bool operator==(const Iterator<OtherType> & rhs) const
    {
        return m_current == rhs.m_current && m_pos == rhs.m_pos;
    }
    template <bool OtherType>
    bool operator!=(const Iterator<OtherType> & rhs) const
    {
        return m_current != rhs.m_current || m_pos != rhs.m_pos;
    }

    reference operator*() const
    {
        return m_current->data[m_pos];
    }

    pointer operator->() const
    {
        return &m_current->data[m_pos];
    }

    Iterator & operator++()
    {
        ++m_pos;
        if (m_pos >= m_current->size() && m_current->right != nullptr) {
            m_pos = 0;
            m_current = m_current->right;
        }
        return *this;
    }

    Iterator operator++(int)
    {
        auto tmp = *this;
        operator++();
        return tmp;
    }

private:
    friend class BPTree;

    Iterator(Leaf * node)
        : m_current(node)
    {
    }

    Iterator(Leaf * node, const size_type & pos)
        : m_current(node)
        , m_pos(pos)
    {
    }

    Leaf * m_current = nullptr;
    size_type m_pos = 0;
};

template <class Key, class Value, std::size_t BlockSize, class Less>
inline bool BPTree<Key, Value, BlockSize, Less>::less(const Key & lhs, const Key & rhs)
{
    return Less()(lhs, rhs);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline BPTree<Key, Value, BlockSize, Less>::BPTree()
    : m_root(new Leaf)
{
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline BPTree<Key, Value, BlockSize, Less>::BPTree(std::initializer_list<std::pair<Key, Value>> list)
    : m_root(new Leaf)
{
    insert(list);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline BPTree<Key, Value, BlockSize, Less>::BPTree(const BPTree & tree)
    : m_root(copy(tree.m_root))
    , m_size(tree.m_size)
{
    redirect();
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::operator=(const BPTree & tree)
        -> BPTree &
{
    m_root = copy(tree.m_root);
    m_size = tree.m_size;
    redirect();
    return *this;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::copy(Node * node)
        -> Node *
{
    if (node->get_child(0) == nullptr) {
        return new Leaf(*dynamic_cast<Leaf *>(node));
    }
    else {
        Inner * inner = new Inner(*dynamic_cast<Inner *>(node));
        for (size_type i = 0; i <= node->size(); ++i) {
            inner->children.emplace_back(copy(node->get_child(i)));
            inner->children.back()->parent = inner;
        }
        return inner;
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::redirect()
{
    if (m_root != nullptr && m_root->get_child(0) != nullptr) {
        Leaf * leaf = find_min_leaf(m_root);
        if (leaf->parent == nullptr) {
            return;
        }
        Inner * node = leaf->parent;
        for (size_type ind = 0; ind < node->size(); ++ind) {
            leaf = dynamic_cast<Leaf *>(node->get_child(ind));
            leaf->right = dynamic_cast<Leaf *>(node->get_child(ind + 1));
        }
        Inner * next = find_right(node);
        while (next != nullptr) {
            leaf = dynamic_cast<Leaf *>(node->get_child(node->size()));
            leaf->right = dynamic_cast<Leaf *>(next->get_child(0));
            node = next;
            next = find_right(node);
            for (size_type ind = 0; ind < node->size(); ++ind) {
                leaf = dynamic_cast<Leaf *>(node->get_child(ind));
                leaf->right = dynamic_cast<Leaf *>(node->get_child(ind + 1));
            }
        }
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline BPTree<Key, Value, BlockSize, Less>::BPTree(BPTree && tree)
    : m_root(new Leaf)
{
    std::swap(m_root, tree.m_root);
    std::swap(m_size, tree.m_size);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::operator=(BPTree && tree)
        -> BPTree &
{
    m_root = new Leaf;
    std::swap(m_root, tree.m_root);
    std::swap(m_size, tree.m_size);
    return *this;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find_min_leaf(Node * node)
        -> Leaf *
{
    while (node->get_child(0) != nullptr) {
        node = node->get_child(0);
    }
    return dynamic_cast<Leaf *>(node);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find_max_leaf(Node * node)
        -> Leaf *
{
    while (node->get_child(node->size()) != nullptr) {
        node = node->get_child(node->size());
    }
    return dynamic_cast<Leaf *>(node);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::begin()
        -> iterator
{
    return iterator(find_min_leaf(m_root));
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::cbegin() const
        -> const_iterator
{
    return const_iterator(find_min_leaf(m_root));
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::begin() const
        -> const_iterator
{
    return cbegin();
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::end()
        -> iterator
{
    Leaf * leaf = find_max_leaf(m_root);
    return iterator(leaf, leaf->size());
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::cend() const
        -> const_iterator
{
    Leaf * leaf = find_max_leaf(m_root);
    return const_iterator(leaf, leaf->size());
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::end() const
        -> const_iterator
{
    return cend();
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline bool BPTree<Key, Value, BlockSize, Less>::empty() const
{
    return m_size == 0;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::size() const
        -> size_type
{
    return m_size;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::clear()
{
    clear(m_root);
    m_root = new Leaf;
    m_size = 0;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::clear(Node * node)
{
    if (node->get_child(0) != nullptr) {
        Inner * inner = dynamic_cast<Inner *>(node);
        for (auto & child : inner->children) {
            clear(child);
        }
    }
    delete node;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::count(const Key & key) const
        -> size_type
{
    return contains(key) ? 1 : 0;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline bool BPTree<Key, Value, BlockSize, Less>::contains(const Key & key) const
{
    iterator lower = lower_bound(key);
    return (lower != end()) && !less(key, lower->first);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::equal_range(const Key & key)
        -> std::pair<iterator, iterator>
{
    return {lower_bound(key), upper_bound(key)};
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::equal_range(const Key & key) const
        -> std::pair<const_iterator, const_iterator>
{
    return {lower_bound(key), upper_bound(key)};
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::search(const Key & key) const
        -> std::pair<Leaf *, size_type>
{
    Node * node = m_root;
    while (node->get_child(0) != nullptr) {
        Inner * inner = dynamic_cast<Inner *>(node);
        size_type ind = std::upper_bound(inner->keys.begin(), inner->keys.end(), key, less) - inner->keys.begin();
        node = node->get_child(ind);
    }
    Leaf * leaf = dynamic_cast<Leaf *>(node);
    size_type pos = std::lower_bound(leaf->data.begin(), leaf->data.end(), key, [](const auto & lhs, const auto & rhs) { return less(lhs.first, rhs); }) - leaf->data.begin();
    return {leaf, pos};
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::lower_bound(const Key & key)
        -> iterator
{
    auto [leaf, pos] = search(key);
    if (pos >= leaf->size() && leaf->right != nullptr) {
        return iterator(leaf->right);
    }
    return iterator(leaf, pos);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::lower_bound(const Key & key) const
        -> const_iterator
{
    auto [leaf, pos] = search(key);
    if (pos >= leaf->size() && leaf->right != nullptr) {
        return const_iterator(leaf->right);
    }
    return const_iterator(leaf, pos);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::upper_bound(const Key & key)
        -> iterator
{
    iterator lower = lower_bound(key);
    if (lower != end() && !less(key, lower->first)) {
        ++lower;
    }
    return lower;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::upper_bound(const Key & key) const
        -> const_iterator
{
    const_iterator lower = lower_bound(key);
    if (lower != cend() && !less(key, lower->first)) {
        ++lower;
    }
    return lower;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find(const Key & key)
        -> iterator
{
    iterator lower = lower_bound(key);
    if (lower != end() && !less(key, lower->first)) {
        return lower;
    }
    return end();
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find(const Key & key) const
        -> const_iterator
{
    const_iterator lower = lower_bound(key);
    if (lower != cend() && !less(key, lower->first)) {
        return lower;
    }
    return cend();
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline Value & BPTree<Key, Value, BlockSize, Less>::at(const Key & key)
{
    iterator lower = lower_bound(key);
    if (lower != end() && !less(key, lower->first)) {
        return lower->second;
    }
    throw std::out_of_range("");
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline const Value & BPTree<Key, Value, BlockSize, Less>::at(const Key & key) const
{
    const_iterator lower = lower_bound(key);
    if (lower != cend() && !less(key, lower->first)) {
        return lower->second;
    }
    throw std::out_of_range("");
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline Value & BPTree<Key, Value, BlockSize, Less>::operator[](const Key & key)
{
    auto [leaf, pos] = search(key);
    iterator it = iterator(leaf, pos);
    if (pos >= leaf->size() || less(key, leaf->data[pos].first)) {
        it = insert(leaf, pos, key, Value());
        ++m_size;
    }
    return it->second;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::insert(const Key & key, const Value & value)
        -> std::pair<iterator, bool>
{
    auto [leaf, pos] = search(key);
    if (pos >= leaf->size() || less(key, leaf->data[pos].first)) {
        iterator it = insert(leaf, pos, key, value);
        ++m_size;
        return {it, true};
    }
    return {iterator(leaf, pos), false};
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::insert(const Key & key, Value && value)
        -> std::pair<iterator, bool>
{
    auto [leaf, pos] = search(key);
    if (pos >= leaf->size() || less(key, leaf->data[pos].first)) {
        iterator it = insert(leaf, pos, key, std::forward<Value>(value));
        ++m_size;
        return {it, true};
    }
    return {iterator(leaf, pos), false};
}

template <class Key, class Value, std::size_t BlockSize, class Less>
template <class ForwardIt>
inline void BPTree<Key, Value, BlockSize, Less>::insert(ForwardIt begin, ForwardIt end)
{
    for (auto i = begin; i != end; ++i) {
        insert(i->first, i->second);
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::insert(std::initializer_list<value_type> list)
{
    for (const auto & i : list) {
        insert(i.first, i.second);
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
template <class K, class V>
inline auto BPTree<Key, Value, BlockSize, Less>::insert(Leaf * leaf, size_type pos, K && key, V && value)
        -> iterator
{
    leaf->data.insert(leaf->data.begin() + pos, std::make_pair(std::forward<K>(key), std::forward<V>(value)));
    if (leaf->size() == m_order) {
        split(leaf);
        if (pos >= m_order / 2) {
            pos -= m_order / 2;
            leaf = leaf->right;
        }
    }
    return iterator(leaf, pos);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::split(Leaf * left)
{
    Leaf * right = new Leaf;
    std::move(left->data.begin() + m_order / 2, left->data.end(), std::back_inserter(right->data));
    left->data.erase(left->data.begin() + m_order / 2, left->data.end());
    right->right = left->right;
    left->right = right;

    if (left->parent == nullptr) {
        Inner * p = new Inner;
        p->children.push_back(left);
        left->parent = p;
        p->keys.push_back(right->data[0].first);
        p->children.push_back(right);
        right->parent = p;
        m_root = p;
    }
    else {
        Inner * p = dynamic_cast<Inner *>(left->parent);
        right->parent = p;
        size_type ind = std::lower_bound(p->keys.begin(), p->keys.end(), right->data[0].first, less) - p->keys.begin();
        p->keys.insert(p->keys.begin() + ind, right->data[0].first);
        p->children.insert(p->children.begin() + ind + 1, right);
        if (p->size() == m_order) {
            split(p);
        }
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::split(Inner * left)
{
    Inner * right = new Inner;
    std::move(left->keys.begin() + m_order / 2 + 1, left->keys.end(), std::back_inserter(right->keys));
    left->keys.erase(left->keys.begin() + m_order / 2, left->keys.end());
    std::move(left->children.begin() + m_order / 2 + 1, left->children.end(), std::back_inserter(right->children));
    left->children.erase(left->children.begin() + m_order / 2 + 1, left->children.end());
    for (auto & i : right->children) {
        i->parent = right;
    }

    if (left->parent == nullptr) {
        Inner * p = new Inner;
        p->children.push_back(left);
        left->parent = p;
        p->children.push_back(right);
        right->parent = p;
        Leaf * rmin = find_min_leaf(right);
        p->keys.push_back(rmin->data[0].first);
        m_root = p;
    }
    else {
        Inner * p = dynamic_cast<Inner *>(left->parent);
        right->parent = p;
        Leaf * rmin = find_min_leaf(right);
        size_type ind = std::lower_bound(p->keys.begin(), p->keys.end(), rmin->data[0].first, less) - p->keys.begin();
        p->keys.insert(p->keys.begin() + ind, rmin->data[0].first);
        p->children.insert(p->children.begin() + ind + 1, right);
        if (p->size() == m_order) {
            split(p);
        }
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find_left(Leaf * leaf)
        -> Leaf *
{
    if (leaf == nullptr || leaf->size() == 0 || leaf->parent == nullptr) {
        return nullptr;
    }
    const Key & key = leaf->get_key(0);
    Inner * parent = leaf->parent;
    Inner * node = parent;
    size_type height = 0;
    size_type pos = 0;
    while (parent != nullptr && pos == 0) {
        node = parent;
        ++height;
        parent = node->parent;
        pos = std::upper_bound(node->keys.begin(), node->keys.end(), key, less) - node->keys.begin();
    }
    if (parent == nullptr && pos == 0) {
        return nullptr;
    }
    if (height < 2) {
        return dynamic_cast<Leaf *>(node->get_child(pos - 1));
    }
    node = dynamic_cast<Inner *>(node->get_child(pos - 1));
    for (; height > 2; --height) {
        node = dynamic_cast<Inner *>(node->get_child(node->size()));
    }
    return dynamic_cast<Leaf *>(node->get_child(node->size()));
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find_left(Inner * node)
        -> Inner *
{
    if (node == nullptr || node->size() == 0 || node->parent == nullptr) {
        return nullptr;
    }
    const Key & key = node->get_key(0);
    Inner * parent = node->parent;
    size_type height = 0;
    size_type pos = 0;
    while (parent != nullptr && pos == 0) {
        node = parent;
        ++height;
        parent = node->parent;
        pos = std::upper_bound(node->keys.begin(), node->keys.end(), key, less) - node->keys.begin();
    }
    if (parent == nullptr && pos == 0) {
        return nullptr;
    }
    node = dynamic_cast<Inner *>(node->get_child(pos - 1));
    for (; height > 1; --height) {
        node = dynamic_cast<Inner *>(node->get_child(node->size()));
    }
    return node;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find_right(Leaf * node)
        -> Leaf *
{
    return node->right;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::find_right(Inner * node)
        -> Inner *
{
    if (node == nullptr || node->size() == 0 || node->parent == nullptr) {
        return nullptr;
    }
    const Key & key = node->get_key(0);
    Inner * parent = node->parent;
    size_type height = 0;
    size_type pos = node->size();
    while (parent != nullptr && pos >= node->size()) {
        node = parent;
        ++height;
        parent = node->parent;
        pos = std::upper_bound(node->keys.begin(), node->keys.end(), key, less) - node->keys.begin();
    }
    if (parent == nullptr && pos >= node->size()) {
        return nullptr;
    }
    node = dynamic_cast<Inner *>(node->get_child(pos + 1));
    for (; height > 1; --height) {
        node = dynamic_cast<Inner *>(node->get_child(0));
    }
    return node;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::erase(const_iterator begin, const_iterator end)
        -> iterator
{
    if (begin == end) {
        return end;
    }
    std::vector<Key> keys;
    for (auto i = begin; i != end; ++i) {
        keys.push_back(i->first);
    }
    for (const Key & key : keys) {
        erase(key);
    }
    return upper_bound(keys.back());
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::erase(const_iterator it)
        -> iterator
{
    if (it == end()) {
        return end();
    }
    Key key = it->first;
    auto [leaf, pos] = search(key);
    bool res = erase(leaf, pos, key);
    if (res) {
        --m_size;
    }
    return upper_bound(key);
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline auto BPTree<Key, Value, BlockSize, Less>::erase(const Key & key)
        -> size_type
{
    auto [leaf, pos] = search(key);
    bool res = erase(leaf, pos, key);
    if (res) {
        --m_size;
        return 1;
    }
    return 0;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline bool BPTree<Key, Value, BlockSize, Less>::erase(Leaf * node, size_type pos, const Key & key)
{
    if (pos >= node->size() || less(key, node->get_key(pos))) {
        return false;
    }
    Leaf * left = find_left(node);
    Leaf * right = find_right(node);

    Key minimum = key;
    node->data.erase(node->data.begin() + pos);
    if (node->size() > 0 && less(node->get_key(0), minimum)) {
        minimum = node->get_key(0);
    }

    if (node->size() >= m_order / 2) {
        if (pos == 0) {
            update(node, key, node->get_key(0));
        }
        return true;
    }

    if (left != nullptr && left->size() > m_order / 2 && !less(minimum, node->parent->get_key(0))) {
        node->data.insert(node->data.begin(), left->data.back());
        left->data.pop_back();
        update(node, minimum, node->get_key(0));
    }
    else if (right != nullptr && right->size() > m_order / 2) {
        node->data.insert(node->data.end(), right->data.front());
        right->data.erase(right->data.begin());
        if (node->size() == 1) {
            update(node, minimum, node->get_key(0));
        }
        update(right, node->data.back().first, right->get_key(0));
    }
    else if (left != nullptr && !less(minimum, node->parent->get_key(0))) {
        merge(left, node, minimum);
    }
    else if (right != nullptr) {
        merge(node, right, right->get_key(0));
    }

    if (m_root->size() == 0 && m_root->get_child(0) != nullptr) {
        m_root = m_root->get_child(0);
        delete m_root->parent;
        m_root->parent = nullptr;
    }

    return true;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::merge(Leaf * left, Leaf * right, const Key & key)
{
    std::move(right->data.begin(), right->data.end(), std::back_inserter(left->data));
    left->right = right->right;
    Inner * parent = right->parent;
    if (parent != nullptr) {
        size_type pos = std::upper_bound(parent->keys.begin(), parent->keys.end(), key, less) - parent->keys.begin();
        erase(parent, pos, key);
    }
    delete right;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::erase(Inner * node, size_type pos, const Key & key)
{
    Inner * left = find_left(node);
    Inner * right = find_right(node);

    Key minimum = key;
    node->keys.erase(node->keys.begin() + pos - 1);
    node->children.erase(node->children.begin() + pos);
    if (node->size() > 0 && less(node->get_key(0), minimum)) {
        minimum = node->get_key(0);
    }

    if (node->size() >= m_order / 2) {
        if (pos == 0) {
            update(node, key, find_min_leaf(node)->get_key(0));
        }
        return;
    }

    if (left != nullptr && left->size() > m_order / 2 && !less(minimum, node->parent->get_key(0))) {
        node->keys.insert(node->keys.begin(), find_min_leaf(node)->get_key(0));
        left->keys.pop_back();
        node->children.insert(node->children.begin(), left->children.back());
        left->children.back()->parent = node;
        left->children.pop_back();
        update(node, minimum, find_min_leaf(node)->get_key(0));
    }
    else if (right != nullptr && right->size() > m_order / 2) {
        node->keys.insert(node->keys.end(), find_min_leaf(right)->get_key(0));
        right->keys.erase(right->keys.begin());
        node->children.insert(node->children.end(), right->children.front());
        right->children.front()->parent = node;
        right->children.erase(right->children.begin());
        if (node->size() == 1) {
            update(node, minimum, find_min_leaf(node)->get_key(0));
        }
        update(right, node->keys.back(), find_min_leaf(right)->get_key(0));
    }
    else if (left != nullptr && !less(minimum, node->parent->get_key(0))) {
        merge(left, node, minimum);
    }
    else if (right != nullptr) {
        merge(node, right, right->get_key(0));
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::merge(Inner * left, Inner * right, const Key & key)
{
    left->keys.push_back(find_min_leaf(right)->get_key(0));
    std::move(right->keys.begin(), right->keys.end(), std::back_inserter(left->keys));
    for (size_type i = 0; i <= right->size(); ++i) {
        left->children.push_back(right->get_child(i));
        left->children.back()->parent = left;
    }
    Inner * parent = right->parent;
    if (parent != nullptr) {
        size_type pos = std::upper_bound(parent->keys.begin(), parent->keys.end(), key, less) - parent->keys.begin();
        erase(parent, pos, key);
    }
    delete right;
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline void BPTree<Key, Value, BlockSize, Less>::update(Node * node, const Key & prev, const Key & key)
{
    if (node == nullptr) {
        return;
    }
    Inner * parent = node->parent;
    if (parent == nullptr) {
        return;
    }
    while (parent != nullptr) {
        size_type pos = std::upper_bound(parent->keys.begin(), parent->keys.end(), prev, less) - parent->keys.begin();
        if (pos > 0) {
            parent->keys[pos - 1] = key;
            return;
        }
        node = parent;
        parent = node->parent;
    }
}

template <class Key, class Value, std::size_t BlockSize, class Less>
inline BPTree<Key, Value, BlockSize, Less>::~BPTree()
{
    clear(m_root);
}
