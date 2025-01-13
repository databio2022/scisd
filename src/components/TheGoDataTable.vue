<script setup lang="ts">
import DataTable from 'datatables.net-vue3'
import DataTablesCore from 'datatables.net'
import { ref, defineProps, computed, onMounted, watch } from 'vue'

import axios from 'axios'

DataTable.use(DataTablesCore)

const columns = [
  { data: 'ID' },
  { data: 'Description' },
  { data: 'BgRatio' },
  { data: 'GeneRatio' },
  { data: 'pvalue' },
  {
    data: 'p.adjust',
    render: function (data, type, row) {
      // Access 'p.adjust' using bracket notation
      return row['p.adjust']
    }
  }
]

// 定义 props
const props = defineProps<{
  data: any[]
}>()

const options = ref({
  pageLength: 15,
  lengthMenu: [
    [15, 20, 30, 40, 50, -1],
    ['15', '20', '30', '40', '50', 'All']
  ],
  order: [[5, 'asc']] //  设置p.adjust默认升序
})

const layout = {
  title: ''
}
const config = {
  responsive: true
}
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
          <th>ID</th>
          <th>Description</th>
          <th>BgRatio</th>
          <th>GeneRatio</th>
          <th>pvalue</th>
          <th>p.adjust</th>
        </tr>
      </thead>
      <tfoot>
        <tr>
          <th>ID</th>
          <th>Description</th>
          <th>BgRatio</th>
          <th>GeneRatio</th>
          <th>pvalue</th>
          <th>p.adjust</th>
        </tr>
      </tfoot>
    </DataTable>
  </div>
</template>

<style>
@import 'datatables.net-dt';
</style>
