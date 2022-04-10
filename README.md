# B+ дерево

## Потребность в блочных структурах данных
Организация памяти в вычислительных системах сложна и иерархична, с существенным ростом стоимости доступа к каждому следующему уровню.
Поэтому выгодными оказываются блочные схемы обращения к памяти - когда за раз читается относительно большой непрерывный блок данных, даже если это повлечёт
за собой дополнительный вычислительные затраты на его обработку.

Особенно актуальна блочная организация для хранения данных на HDD - скорость случайного доступа у дисков очень низка и выгода от чтения больших блоков
весьма существенна.

Пример:
* Процессор с частотой 1 Гц
* Скорость доступа к кешу процессора - примерно 10 тактов = 10 нс
* Скорость доступа к ОЗУ - примерно 100 тактов = 100 нс
* Скорость доступа к HDD (время перемещения головки) - примерно 10 мс - медленнее ОЗУ в 10^5 раз

## Идея
B+ дерево - это дерево с более высокой степенью ветвления (порядком), нежели двоичное. Внутренние узлы дерева порядка N могут иметь от ceil(N/2) до N потомков.
Порядок дерева - не меньше 3.
Корень дерева может иметь от 2 до N потомков.

Внутренние узлы дерева хранят только ключи и ссылки на потомков.

Листья дерева содержат ключи и значения, каждый лист может содержать от N/2 до N пар ключ-значение. Листья связаны в односвязный список, что позволяет эффективно
реализовывать поиск интервалов.

Подробнее можно ознакомиться с идеей B+ дерева в [статье википедии](https://en.wikipedia.org/wiki/B%2B_tree).
Другое описание см. [здесь](http://www.cburch.com/cs/340/reading/btree/index.html).

Основных отличия от двоичных деревьев несколько:
* узлы гораздо большего размера, хранят много ключей
* внутренние узлы хранят только ключи и ссылки на потомков, любой поиск сопряжён со спуском до листа
* листья связаны в односвязный список, что позволяет эффективно получать интервалы

Обычно размер узла делают достаточно существенным, чтобы читать данные большими блоками. Например, если на хранение ключа и ссылки на потомка требуется 32 байта,
размер узла - 4096 байт, то порядок дерева = 128. Если требуется хранить в таком дереве 10^7 элементов, то высота его будет варьироваться от 4 до 5 - то есть,
любой поиск будет сопряжён максимум с 5 чтениями. Сравните с бинарным деревом - у идеально сбалансированного дерева на 10^7 элементов высота будет равна 23.

## Постановка задачи
Требуется реализовать шаблонный ассоциативный массив с интерфейсом, во многом подобным `std::map`.

Предполагается, что все хранимые ключи уникальны (как и с `std::map`).

Нужно реализовать основные функции типичной коллекции:
* проверка размера и очистка
* провера наличия ключа
* поиск элемента по ключу
* прямой доступ к элементу по ключу
* добавление и удаление элементов
* поддержка задания произвольного компаратора

В основе реализации должно лежать B+ дерево. Его порядок не фиксирован и выбирается исходя из размера хранимых типов так, чтобы обеспечить фиксированный
размер блока, задаваемый шаблонным параметром класса дерева, значение по умолчанию - 4096 байт.
Предельная заполненной блока должна быть максимально возможной, например, если для обслуживания 1 ключа во внутреннем блоке
требуется 100 байт, то порядок дерева будет равен 40.

Допускается во внутренних узлах иметь дополнительные накладные расходы в виде размер ключа + размер указателя (т.е. при эффективной вместимости узла в `K` ключей,
отводить место под `K+1`), это может упростить реализацию.
Помимо этого предполагается, что внутренние узлы имеют место для хранения `K` ключей и `K+1` указателей, а также хранят количество присутствующих ключей.
Таким образом, `K` можно определить как `(BlockSize - sizeof(std::size_t) - sizeof(void *)) / (sizeof(Key) + sizeof(void *))`.
Листья имеют тот же размер, который используют для хранения пар ключ-значение, количества присутствующих пар и указателя на следующий лист.

Возможность задания произвольных типов ключей и значений нарушает идею "блочной" структуры данных: к примеру, `std::string` выделяет память под своё значение
самостоятельно и на это сложно повлиять. Тем не менее, по сравнению с `std::map`, такое B+ дерево всё-таки будет "блочным" - доступ к одному узлу даст множество
ключей, а не один.

Итераторы коллекции должны иметь категорию forward.

Допустимы следующие упрощения:
* можно игнорировать служебные расходы на организацию классов (например, для поддержки виртуальных функций); достаточно, если в заданный размер блока укладывается
  всё явно определённое содержимое узла
* можно игнорировать то, что некоторые классы ключей или значений могут самостоятельно выделять память; достаточно упаковать в заданный размер блока сами объекты
  ключей и значений (и иные вспомогательные поля классов узлов)
* разное количество ключей во внутренних узлах и листьях: размер блока распространяется на оба типа узлов, но если внутренние узлы помимо ключей содержат ссылки на
  потомков, то листья содержат помимо ключей - объекты значений, размер которых может отличаться от затрат на представление ссылок

## Требуемый интерфейс
```c++
template <class Key, class Value, std::size_t BlockSize, class Less = std::less<Key>>
class BPTree
{
public:
    using key_type = Key;
    using mapped_type = Value;
    using value_type = std::pair<Key, Value>; // NB: a digression from std::map
    using reference = value_type &;
    using const_reference = const value_type &;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using size_type = std::size_t;

    using iterator = ...;
    using const_iterator = ...;

    BPTree(std::initializer_list<std::pair<Key, Value>>);
    BPTree(const BPTree &);
    BPTree(BPTree &&);

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
    Value & operator[] (const Key &);

    std::pair<iterator, bool> insert(const Key &, const Value &); // NB: a digression from std::map
    std::pair<iterator, bool> insert(const Key &, Value &&); // NB: a digression from std::map
    template <class ForwardIt>
    void insert(ForwardIt begin, ForwardIt end);
    void insert(std::initializer_list<value_type>);
    iterator erase(const_iterator);
    iterator erase(const_iterator, const_iterator);
    size_type erase(const Key &);
};
```

## Рекомендации
Данная структура данных весьма проста идейно, но её реализация оказывается довольно сложной. Причём эта сложность не алгоритмическая, а во многом связана с организацией
структуры на конкретном языке. Самая сложная операция - удаления. Попробуйте для начала отключить тесты, требующие операций удаления и реализовать только вставку.

У некоторых авторов возникает желание добавить горизонтальные связи между внутренними узлами дерева - лучше от этого воздержаться, т.к. только операция удаления может
получить от этого какой-то выигрыш, но для удаления вам всё равно понадобятся сложные вертикальные манипуляции и соседа можно получить через родителя.

Если допустить возможность временного переполнения внутренних узлов при вставке (см. выше про допустимые накладные расходы), то операцию вставки во внутренний узел
можно сильно упростить выделив две раздельные стадии - вставку и расщепление.

Говорят, что на викиконспектах ИТМО в описании операций над B+ деревом есть ошибки или неточности.

При возникновении проблем, попробуйте воспроизвести каждую из них по отдельности, написав минимальный юнит-тест (за хорошие тесты вы можете получить дополнительные баллы).
