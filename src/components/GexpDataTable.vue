<script setup lang="ts">
import DataTable from 'datatables.net-vue3'
import DataTablesCore from 'datatables.net'
import { ref, defineProps } from 'vue'

DataTable.use(DataTablesCore)

// 定义 props
const props = defineProps<{
  data: any[]
  gene: String
}>()

const columns = [
  { data: 'gene' },
  {
    data: null,
    title: 'targetgene',
    render: function (data, type, row, meta) {
      return props.gene
    }
  },
  {
    data: 'cor'
  },
  {
    data: 'pvalue'
  }
]

const options = ref({
  pageLength: 15,
  lengthMenu: [
    [15, 20, 30, 40, 50, -1],
    ['15', '20', '30', '40', '50', 'All']
  ],
  ordering: true,
  order: [[2, 'desc']]
})
</script>

<template>
  <div>
    <DataTable
      :columns="columns"
      :data="props.data"
      class="display"
      :options="options"
      width="100%"
    >
      <thead>
        <tr>
          <th>corrgene</th>
          <th>targetgene</th>
          <th>Pearson R</th>
          <th>P Value</th>
        </tr>
      </thead>
      <tfoot>
        <tr>
          <th>corrgene</th>
          <th>targetgene</th>
          <th>Pearson R</th>
          <th>P Value</th>
        </tr>
      </tfoot>
    </DataTable>
  </div>
</template>

<style>
@import 'datatables.net-dt';
</style>
