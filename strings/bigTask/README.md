Реализации алгоритмов [Фараха](http://www.cs.rutgers.edu/~farach/pubs/Suffix.pdf) и [индуцированной сортировки](https://local.ugene.unipro.ru/tracker/secure/attachment/12144/Linear+Suffix+Array+Construction+by+Almost+Pure+Induced-Sorting.pdf)

В файлах `<algorithm name>/toSubmit<online judge name>.cpp` лежат коды посылок в соответствующую систему (получены конкатенацией всех прочих файлов, удалением лишних `include`-ов и дописыванием реализации требуемой задачи с использованием указанного алгоритма)

Задачи:
* [1706. Шифровка 2](http://acm.timus.ru/problem.aspx?space=1&num=1706) на Тимусе
* [111789. Различные подстроки](http://informatics.mccme.ru/moodle/mod/statements/view3.php?chapterid=111789) на информатиксе


UPD:

Если скомпилировать `inducedSorting.cpp` с `-O2` (построить суфмас на `1e7` элементов), то будет работать `4.406181` секунд. Успех.
