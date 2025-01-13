<template>
  <div class="container mx-auto mt-4">
    <div class="card bg-base-100 border">
      <div class="p-4 bg-gray-200">
        <h2 class="text-lg font-bold">Find Conserved Marker Genes Across Datasets</h2>
      </div>
      <div class="card-body">
        <div class="flex space-x-4 justify-center">
          <div class="relative">
            <label class="input input-bordered flex items-center gap-2 max-w-sm">
              <input
                type="text"
                class="grow"
                placeholder="Search gene, e.g. ACTB"
                v-model="selectedGene"
              />
              <svg
                xmlns="http://www.w3.org/2000/svg"
                viewBox="0 0 16 16"
                fill="currentColor"
                class="h-4 w-4 opacity-70"
              >
                <path
                  fill-rule="evenodd"
                  d="M9.965 11.026a5 5 0 1 1 1.06-1.06l2.755 2.754a.75.75 0 1 1-1.06 1.06l-2.755-2.754ZM10.5 7a3.5 3.5 0 1 1-7 0 3.5 3.5 0 0 1 7 0Z"
                  clip-rule="evenodd"
                />
              </svg>
            </label>
            <div
              v-if="!isGeneSelected && filteredGenes.length"
              class="absolute top-full left-0 w-full z-10 bg-white border shadow-md"
            >
              <VirtualList :list="filteredGenes" @gene-selected="handleGeneSelected" />
            </div>
          </div>
          <button class="btn btn-accent" @click="getbarplot">Search</button>
        </div>
        <!-- waiting animation-->
        <div class="flex items-center justify-center m-24 p-24" v-if="isLoading">
          <div
            class="animate-spin rounded-full h-12 w-12 border-t-4 border-blue-500 border-opacity-75"
          ></div>
        </div>
        <div v-if="loading">
          <BarPlotComponent :data="bardata" :layout="layout"></BarPlotComponent>
        </div>
      </div>
    </div>

    <div class="card bg-base-100 border my-4">
      <div class="p-4 bg-gray-200">
        <h2 class="text-lg font-bold">Find Conserved Disease-related Genes Across Datasets</h2>
      </div>
      <div class="card-body">
        <div>
          <div class="flex space-x-4 justify-between">
            <div class="flex space-x-4">
              <select class="select select-bordered w-full max-w-xs" v-model="selectedCell">
                <option v-for="type in celltypes" :key="type" :value="type">
                  {{ type }}
                </option>
              </select>
              <div class="relative">
                <label class="input input-bordered flex items-center gap-2 max-w-sm">
                  <input
                    type="text"
                    class="grow"
                    placeholder="Search gene, e.g. ACTB"
                    v-model="XselectedGene"
                  />
                  <svg
                    xmlns="http://www.w3.org/2000/svg"
                    viewBox="0 0 16 16"
                    fill="currentColor"
                    class="h-4 w-4 opacity-70"
                  >
                    <path
                      fill-rule="evenodd"
                      d="M9.965 11.026a5 5 0 1 1 1.06-1.06l2.755 2.754a.75.75 0 1 1-1.06 1.06l-2.755-2.754ZM10.5 7a3.5 3.5 0 1 1-7 0 3.5 3.5 0 0 1 7 0Z"
                      clip-rule="evenodd"
                    />
                  </svg>
                </label>
                <div
                  v-if="!XisGeneSelected && XfilteredGenes.length"
                  class="absolute top-full left-0 w-full z-10 bg-white border shadow-md"
                >
                  <VirtualList :list="XfilteredGenes" @gene-selected="XhandleGeneSelected" />
                </div>
              </div>
              <button class="btn btn-accent" @click="gettable">Search</button>
            </div>
            <div>
              <!-- Open the modal using ID.showModal() method -->
              <button class="btn btn-accent" @click="showModalAndPlot">Plot</button>
              <dialog id="my_modal_1" class="modal" ref="myModal">
                <div class="modal-box max-w-6xl">
                  <h3 class="text-lg font-bold">
                    The Figure show fold changes of the gene of interest across datasets,red
                    represents up-regualtion in disease relative to normal, while blue represents
                    down-regulation.
                  </h3>
                  <div v-if="showtable">
                    <BarPlotComponent :data="xbardata" :layout="xlayout"></BarPlotComponent>
                  </div>
                  <div class="modal-action">
                    <form method="dialog">
                      <!-- if there is a button in form, it will close the modal -->
                      <button class="btn" @click="fclose">Close</button>
                    </form>
                  </div>
                </div>
              </dialog>
            </div>
          </div>
          <!-- waiting animation-->
          <div class="flex items-center justify-center m-24 p-24" v-if="isLoading">
            <div
              class="animate-spin rounded-full h-12 w-12 border-t-4 border-blue-500 border-opacity-75"
            ></div>
          </div>

          <div v-if="loading">
            <TheCdegDataTable :data="tabledata"></TheCdegDataTable>
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import TheCdegDataTable from '../components/TheCdegDataTable.vue'
import BarPlotComponent from '../components/BarplotChart.vue'
import axios from 'axios'
import VirtualList from '../components/VirtualList.vue'

export default {
  components: {
    BarPlotComponent,
    TheCdegDataTable,
    VirtualList
  },
  data() {
    return {
      loading: false,
      isLoading: true,

      selectedGene: 'ACTB',

      tabledata: [],

      bardata: [],
      layout: null,
      selectedCell: 'Basal keratinocytes',
      celltypes: [],

      isGeneSelected: false,
      XselectedGene: 'ACTB',
      XisGeneSelected: false,
      first: false,

      genes: [],
      degenes: [],
      xlayout: null,
      xbardata: [],
      showtable: false
    }
  },
  mounted() {
    this.fetchData()
  },

  computed: {
    filteredGenes() {
      if (!this.selectedGene || this.selectedGene === 'ACTB') {
        // 或者 if (this.selectedGene.trim() === '')
        return [] // 输入框为空时返回空数组
      }
      return this.genes.filter((gene) =>
        gene.toLowerCase().includes(this.selectedGene.toLowerCase())
      )
    },
    XfilteredGenes() {
      if (!this.XselectedGene || this.XselectedGene === 'ACTB') {
        // 或者 if (this.selectedGene.trim() === '')
        return [] // 输入框为空时返回空数组
      }
      return this.degenes.filter((gene) =>
        gene.toLowerCase().includes(this.XselectedGene.toLowerCase())
      )
    }
  },

  watch: {
    selectedGene(newValue) {
      if (newValue === '') {
        this.isGeneSelected = false // 输入框为空时显示列表
      }
    },
    XselectedGene(newValue) {
      if (newValue === '') {
        this.XisGeneSelected = false // 输入框为空时显示列表
      }
    }
  },

  methods: {
    handleGeneSelected(gene) {
      this.selectedGene = gene
      this.isGeneSelected = true
    },
    XhandleGeneSelected(gene) {
      this.XselectedGene = gene
      this.XisGeneSelected = true
    },
    showModalAndPlot() {
      this.showtable = true
      this.$refs.myModal.showModal()
    },
    fclose() {
      this.showtable = false
    },
    gettable() {
      //table data
      axios
        .post('http://scisd.databio1.com/api/deg_conserve', {
          gene: this.XselectedGene,
          celltype: this.selectedCell
        })
        .then((response1) => {
          this.tabledata = response1.data
          //barplot

          const outputData = this.tabledata.reduce(
            (acc, curr) => {
              acc['x'].push(curr.final_vs)
              acc['y'].push(Math.abs(curr.avg_log2FC))
              acc['marker']['color'].push(
                curr.avg_log2FC > 0 ? 'rgba(219, 64, 82, 0.7)' : 'rgba(54, 162, 235, 0.7)'
              )
              return acc
            },
            { x: [], y: [], marker: { color: [] } }
          )
          outputData['type'] = 'bar'

          this.xbardata = [outputData]

          const layout = {
            title: 'Consered disease-related genes across datasets: ' + this.XselectedGene,

            yaxis: { title: 'avg_log2FC' },
            xaxis: {
              tickangle: -45, //  关键代码：设置倾斜角度为-45度
              automargin: true
            }
          }
          this.xlayout = layout

          let data = [outputData]
          // 提取 x 和 y 值
          const x = data[0].x
          const y = data[0].y

          // 创建一个包含 x 和 y 值对的数组
          const combined = x.map((xValue, index) => ({ x: xValue, y: y[index] }))

          // 按 y 值降序排序
          combined.sort((a, b) => b.y - a.y)

          // 从排序后的数组中提取 x 和 y 值
          data[0].x = combined.map((pair) => pair.x)
          data[0].y = combined.map((pair) => pair.y)
        })
        .catch((error) => {
          console.error(error)
        })
    },
    getbarplot() {
      this.isLoading = true
      this.loading = false
      // 使用 axios 发送 POST 请求
      axios
        .post('http://scisd.databio1.com/api/marker_conserve', {
          gene: this.selectedGene
        })
        .then((response) => {
          this.substatus = false
          // 处理服务器返回的数据
          //console.log(response.data)

          const outputData = response.data.reduce(
            (acc, curr) => {
              acc['x'].push(curr.celltype)
              acc['y'].push(curr.frequency)
              return acc
            },
            { x: [], y: [] }
          )
          outputData['type'] = 'bar'
          outputData['marker'] = {
            color: 'rgba(219, 64, 82, 0.7)',
            line: {
              color: 'rgba(219, 64, 82, 1.0)',
              width: 1
            }
          }

          this.bardata = [outputData]
          const layout = {
            title: 'Consered marker genes across datasets: ' + this.selectedGene,

            yaxis: { title: 'Frequency' }
          }
          this.layout = layout

          let data = [outputData]
          // 提取 x 和 y 值
          const x = data[0].x
          const y = data[0].y

          // 创建一个包含 x 和 y 值对的数组
          const combined = x.map((xValue, index) => ({ x: xValue, y: y[index] }))

          // 按 y 值降序排序
          combined.sort((a, b) => b.y - a.y)

          // 从排序后的数组中提取 x 和 y 值
          data[0].x = combined.map((pair) => pair.x)
          data[0].y = combined.map((pair) => pair.y)

          this.isLoading = false
          this.loading = true
        })
        .catch((error) => {
          console.error(error)
        })
    },

    fetchData() {
      if (this.selectedGene) {
        // 使用 axios 发送 POST 请求
        axios
          .post('http://scisd.databio1.com/api/marker_conserve', {
            gene: this.selectedGene
          })
          .then((response) => {
            this.substatus = false
            // 处理服务器返回的数据
            //console.log(response.data)

            const outputData = response.data.reduce(
              (acc, curr) => {
                acc['x'].push(curr.celltype)
                acc['y'].push(curr.frequency)
                return acc
              },
              { x: [], y: [] }
            )
            outputData['type'] = 'bar'
            outputData['marker'] = {
              color: 'rgba(219, 64, 82, 0.7)',
              line: {
                color: 'rgba(219, 64, 82, 1.0)',
                width: 1
              }
            }
            this.bardata = [outputData]
            const layout = {
              title: 'Consered marker genes across datasets: ' + this.selectedGene,

              yaxis: { title: 'Frequency' }
            }
            this.layout = layout

            let data = [outputData]
            // 提取 x 和 y 值
            const x = data[0].x
            const y = data[0].y

            // 创建一个包含 x 和 y 值对的数组
            const combined = x.map((xValue, index) => ({ x: xValue, y: y[index] }))

            // 按 y 值降序排序
            combined.sort((a, b) => b.y - a.y)

            // 从排序后的数组中提取 x 和 y 值
            data[0].x = combined.map((pair) => pair.x)
            data[0].y = combined.map((pair) => pair.y)

            //get genelist
            axios
              .get('http://scisd.databio1.com/api/marker_conserve_gene')
              .then((responseb) => {
                const x = JSON.parse(responseb.data)
                const geneArray1 = x.map((item) => Object.values(item)[0])
                this.genes = geneArray1
                //console.log(geneArray1)
                axios
                  .get('http://scisd.databio1.com/api/deg_conserve_gene')
                  .then((response1) => {
                    const x = JSON.parse(response1.data)
                    const geneArray1 = x.map((item) => Object.values(item)[0])
                    this.degenes = geneArray1
                    //console.log(geneArray1)
                    //celltype
                    axios
                      .get('http://scisd.databio1.com/api/de_conserve_celltype')
                      .then((response1) => {
                        const x = JSON.parse(response1.data)
                        const geneArray1 = x.map((item) => Object.values(item)[0])
                        this.celltype = geneArray1[0]
                        this.celltypes = geneArray1

                        //console.log(geneArray1)
                        //table data
                        axios
                          .post('http://scisd.databio1.com/api/deg_conserve', {
                            gene: this.selectedGene,
                            celltype: this.selectedCell
                          })
                          .then((response1) => {
                            //console.log(response1.data)
                            this.tabledata = response1.data
                            const outputData = this.tabledata.reduce(
                              (acc, curr) => {
                                acc['x'].push(curr.final_vs + ',' + curr.dataset)
                                acc['y'].push(curr.avg_log2FC)
                                return acc
                              },
                              { x: [], y: [] }
                            )
                            outputData['type'] = 'bar'
                            outputData['marker'] = {
                              color: 'rgba(219, 64, 82, 0.7)',
                              line: {
                                color: 'rgba(219, 64, 82, 1.0)',
                                width: 1
                              }
                            }

                            this.xbardata = [outputData]
                            const layout = {
                              height: 600,
                              title:
                                'Consered disease-related genes across datasets: ' +
                                this.XselectedGene,

                              yaxis: { title: 'avg_log2FC' },
                              xaxis: {
                                tickangle: -45, //  关键代码：设置倾斜角度为-45度
                                automargin: true
                              }
                            }
                            this.xlayout = layout

                            let data = [outputData]
                            // 提取 x 和 y 值
                            const x = data[0].x
                            const y = data[0].y

                            // 创建一个包含 x 和 y 值对的数组
                            const combined = x.map((xValue, index) => ({ x: xValue, y: y[index] }))

                            // 按 y 值降序排序
                            combined.sort((a, b) => b.y - a.y)

                            // 从排序后的数组中提取 x 和 y 值
                            data[0].x = combined.map((pair) => pair.x)
                            data[0].y = combined.map((pair) => pair.y)
                          })
                          .catch((error) => {
                            console.error(error)
                          })
                      })
                      .catch((error) => {
                        console.error(error)
                      })
                  })
                  .catch((error) => {
                    console.error(error)
                  })
              })
              .catch((error) => {
                console.error(error)
              })

            this.isLoading = false
            this.loading = true
          })
          .catch((error) => {
            console.error(error)
          })
      }
    }
  }
}
</script>
