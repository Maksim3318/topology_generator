# Topology Generator

Утилита для генерации и анализа топологий сетей на основе различных алгоритмов.

## Требования

* C++ компилятор (с поддержкой C++17 и выше)
* CMake (версия 3.0 или выше)

## Сборка

Сборка в режиме Release осуществляется через CMake:

```bash
mkdir -p build
cmake -S . -B build -DCMAKE_RUNTIME_OUTPUT_DIRECTORY="$(pwd)"
cmake --build build --config Release
```

После успешной сборки исполняемый файл `topology_generator` появится в директории основной директории проекта.

## Запуск и режимы работы

Утилита поддерживает два режима работы:

* **params** — вычисление метрик топологии
* **net** — генерация конфигурации для Anynet

**Формат запуска:**

```bash
./topology_generator <mode> [filename] [delays] <topology> [parameters...]
```

* `<mode>` — один из `params` или `net`
* `<topology>` — название топологии (чувствительно к регистру)
* `[parameters...]` — список параметров, зависящий от выбранной топологии (см. таблицу ниже)

## Формат входных данных и конфигурация

* Конфигурация приложения загружается из файла `config.conf` в корне репозитория.
* Для режима `net` дополнительно указывается имя выходного файла с сетевой конфигурацией и задержки `router_to_node`, `node_to_router`, `router_to_router`:

  ```bash
  ./topology_generator net <output_filename> <delay r_n n_r r_r>  <topology> [parameters...]
  ```

## Таблица топологий и их параметров

| Topology | Parameters |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| **3DDeBruijn**              | `d (int)` – размерность De Bruijn графа; `k (int)` – размер алфавита; `layers (int)` – число слоёв |
| **3DMesh**                  | `x_dim (int)` – размер по оси X; `y_dim (int)` – размер по оси Y; `z_dim (int)` – размер по оси Z |
| **3DTorus**                 | `x_dim (int)` – размер по оси X; `y_dim (int)` – размер по оси Y; `z_dim (int)` – размер по оси Z |
| **BFT**                     | `radix (int)` – количество детей на каждом уровне; `levels (int)` – число уровней|
| **Banyan**                  | `num_cores (int)` – количество ядер|
| **Benes**                   | `num_cores (int)` – количество ядер|
| **Butterfly**               | `num_cores (int)` – количество ядер|
| **C2Mesh**                  | `rows (int)` – число строк; `cols (int)` – число столбцов |
| **C2Torus**                 | `rows (int)` – число строк; `cols (int)` – число столбцов |
| **CBPMesh**                 | `rows (int)` – число строк; `cols (int)` – число столбцов |
| **CBPTorus**                | `rows (int)` – число строк; `cols (int)` – число столбцов |
| **CCN\_HIN**                | `sizes (list<int>)` – список размеров для каждого уровня|
| **Circulant**               | `num_vertices (int)` – число вершин; `generator_offsets... (int...)` – набор образующих (список чисел через пробел) |
| **Clos**                    | `num_terminals (int)` – количество ядер; `radix (int)` – степень узлов центрального слоя|
| **Cnoc**                    | `groups (int)` – число групп; `terminals_per_group (int)` – количество ядер в группе; `radix (int)` – степень узлов центрального слоя    |
| **CubeConnectedCirculants** | `hyperDim (int)` – размерность гиперкуба; `circBase (int)` – база мультипликативного циркулянта; `circExp (int)` – степень мультипликативного циркулянта |
| **CubeConnectedCycles**     | `hyperDim (int)` – размерность гиперкуба |
| **CubeTreeHybrid**          | `levels (int)` – количество уровней  дерева |
| **DCM**                     | `rows (int)` – число строк; `cols (int)` – число столбцов |
| **DIMB**                    | `d (int)` – размерность De Bruijn графа; `k (int)` – размер алфавита |
| **DMesh**                   | `rows (int)` – число строк; `cols (int)` – число столбцов|
| **DTorus**                  | `rows (int)` – число строк; `cols (int)` – число столбцов|
| **DeBruijn**                | `d (int)` – размерность De Bruijn графа; `k (int)` – размер алфавита |
| **Delta**                   | `num_cores (int)` – количество ядер |
| **DiagonalToroidalMesh**    | `k1 (int)` – смещение по главной диагонали; `k2 (int)` – смещение по побочной диагонали |
| **DiametricalMesh**         | `rows (int)` – число строк; `cols (int)` – число столбцов|
| **Dragonfly**               | `p (int)` – число роутеров в группе; `a (int)` – число групп; `g (int)` – глобальных связей |
| **FBFT**                    | `radix (int)` – степень вершин; `levels (int)` – число уровней |
| **FT**                      | `radix (int)` – степень вершин; `levels (int)` – число уровней |
| **FibonacciCube**           | `n (int)` – размерность гиперкуба |
| **FlattenedButterfly**      | `radix (int)` – степень вершин; `layers (int)` – число слоёв |
| **Flip**                    | `num_cores (int)` – количество ядер |
| **GFT**                     | `h (int)` – высота; `m (int)` – ширина; `w (int)` – глубина |
| **Gaussian**                | `a (int)`; `b (int)`; `N = a^2 + b^2`  |
| **H3DMesh**                 | `m (int)` – число строк; `n (int)` – число столбцов; `L (int)` – число слоёв |
| **H3DTorus**                | `m (int)` – строки; `n (int)` – столбцы; `L (int)` – слоёв|
| **HERT**                    | `k (int)` – количество колец; `m (int)` – количество групп в кольце |
| **HSMBFT**                  | `radix (int)` – основание; `depth (int)` – глубина|
| **HTN**                     | `m (int)` – строки; `n (int)` – столбцы; `L (int)` – слоёв |
| **HexStarMesh**             | `rows (int)` – строки; `cols (int)` – столбцы; `d (int)` – глубина звезды (1 или 2) |
| **HierarchicalClique**      | `radix (int)` – степень вершин; `levels (int)` – уровни |
| **HoneycombMesh**           | `t (int)` – `N = 6t^2` |
| **HoneycombTorus**          | `t (int)` – `N = 6t^2` |
| **HyperCirculant**          | `num_vertices (int)` – число вершин; `a_1,a_2...a_n;b_1,b_2...b_n;...;x_1,x_2...x_n (int...)` – образующие |
| **HyperX**                  | `dimensions (list<int>)` – размеры по каждому измерению |
| **Hypercube**               | `dimensions (int)` – число измерений |
| **Hypermesh**               | `dimensions (list<int>)` – размеры по измерениям |
| **Hypertorus**              | `dimensions (list<int>)` – размеры по измерениям|
| **KAryNCubes**              | `size_per_dim (int)` – размерность по измерению; `num_dimensions (int)` – число измерений |
| **L2Star**                  | `rows (int)` – число строк; `cols (int)` – число столбцов |
| **LNetwork**                | `a (int)` – A; `b (int)` – B; `c (int)` – C; `d (int)` – D |
| **MCC**                     | `rows (int)` – число строк; `cols (int)` – число столбцов |
| **MH3DT**                   | `m (int)` – строки; `n (int)` – столбцы; `L (int)` – слоёв; `q (int)` – ветвлений |
| **MTESH**                   | `sizes (list<int>)` – размеры каждого измерения |
| **MTN**                     | `m (int)` – размер; `level (int)` – уровень; `q (int)` – ветвлений |
| **Mesh**                    | `rows (int)` – число строк; `cols (int)` – число столбцов |
| **MeshMesh**                | `meshSide (int)` – размер mesh |
| **MeshOfSpidergon**         | `n (int)` – число вершин spidergon; `meshSide (int)` – размер стороны |
| **MeshOfStars**             | `leavesPerHub (int)` – листьев на хаб; `meshSide (int)` – размер стороны |
| **MeshOfTrees**             | `mesh_k (int)` – размер дерева |
| **MeshRCTM**                | `meshSide (int)` – размер стороны |
| **Metacube**                | `k (int)` – размерность; `m (int)` – глубина |
| **Midimew**                 | `N (int)` – количество вершин; `b (int)` – параметр циркулянта |
| **MultiplicativeCirculant** | `base (int)` – основание; `degree (int)` – степень |
| **OCBPMesh**                | `rows (int)` – число строк; `cols (int)` – число столбцов |
| **OCBPTorus**               | `rows (int)` – число строк; `cols (int)` – число столбцов |
| **Omega**                   | `num_cores (int)` – число ядер |
| **Paley**                   | `num_vertices (int)` – число вершин; `degree (int)` – степень|
| **Pyramid**                 | `radix (int)` – степень вершин дерева; `levels (int)` – уровни |
| **RCR**                     | `k (int)` – размерность; `r (int)` – радиус; `j (int)` – смещение                                                   |
| **RDT**                     | `rows (int)` – строки; `cols (int)` – столбцы; `R (int)` – радиус; `n (int)` – порядок |
| **RNT**                     | `k (int)` – размерность; `L (int)` – слоёв; `layers (int)` – количество слоёв                                       |
| **Ricobit**                 | `K (int)` – количество колец |
| **Ring**                    | `num_vertices (int)` – количество вершин |
| **RingMesh**                | `ringSize (int)` – размер кольца |
| **RingRCTM**                | `ringSize (int)` – размер кольца |
| **SDMesh**                  | `rows (int)` – число строк; `cols (int)` – число столбцов |
| **SDTorus**                 | `rows (int)` – число строк; `cols (int)` – число столбцов |
| **SEM**                     | `R (int)` – размер стороны mesh |
| **SMBFT**                   | `radix (int)` – степень вершин; `depth (int)` – глубина |
| **SMesh**                   | `x_dim (int)` – размер X; `y_dim (int)` – Y; `z_dim (int)` – Z |
| **Spidergon**               | `num_vertices (int)` – количество вершин |
| **StarMesh**                | `rows (int)` – число строк; `cols (int)` – число столбцов |
| **TESH**                    | `sizes (list<int>)` – размеры (числа через пробел) |
| **THIN**                    | `L (int)` – количество уровней |
| **TMesh**                   | `rows (int)` – число строк; `cols (int)` – число столбцов |
| **Torus**                   | `rows (int)` – число строк; `cols (int)` – число столбцов |
| **TwistedTorus**            | `width (int)` – ширина; `height (int)` – высота; `twist (int)` – степень скручивания |
| **WKRecursive**             | `k (int)` – размерность; `L (int)` – количество уровней |
| **XDMesh**                  | `rows (int)` – число строк; `cols (int)` – число столбцов|
| **XDTorus**                 | `rows (int)` – число строк; `cols (int)` – число столбцов|
| **XNetwork**                | `rows (int)` – число строк; `cols (int)` – число столбцов|
| **XTree**                   | `radix (int)` – степень вершни; `levels (int)` – уровни |
| **ZMesh**                   | `rows (int)` – число строк; `cols (int)` – число столбцов |

## Примеры

* Вычислить параметры 4-мерного гиперквадрата:

  ```bash
  ./topology_generator params Hypercube 4
  ```

* Сгенерировать конфигурацию для топологии `Circulant` с 1024 вершинами и смещениями 1, 144, 258, 276:

  ```bash
  ./topology_generator net out.net 1 2 3 Circulant 1024 1 144 258 276
  ```

---

Конфигурационные параметры (`config.conf`) и дополнительные режимы управления поведением утилиты см. в соответствующем файле.
